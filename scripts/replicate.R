source('scripts/common.R')
x <- load_analysis('discovery'); r <- load_analysis('replication')
a <- read_tsv('results/summary/signal_associations.tsv')
candidate <- a[a$scan_significant, ]; candidate <- candidate[order(candidate$p), ]
Z <- ancestry_center(x$G, x$pheno$ancestry_group)
leads <- character(); membership <- data.frame(lead = character(), variant = character(), r2 = numeric())
remaining <- candidate$variant
while (length(remaining)) {
  lead <- remaining[1]; leads <- c(leads, lead)
  j <- match(lead, x$variants$variant)
  near <- which(x$variants$chr == x$variants$chr[j] & abs(x$variants$pos - x$variants$pos[j]) <= cfg$clump_window_bp)
  r2 <- as.vector(cor(Z[, j], Z[, near, drop = FALSE])^2)
  linked <- near[is.finite(r2) & r2 >= cfg$clump_r2]
  membership <- rbind(membership, data.frame(lead, variant = x$variants$variant[linked], r2 = r2[is.finite(r2) & r2 >= cfg$clump_r2]))
  remaining <- setdiff(remaining, x$variants$variant[linked])
}
write_tsv(membership, 'results/replication/clumps.tsv')
rows <- list(); diagnostics <- list(); site_rows <- list()
for (lead in leads) {
  drow <- a[match(lead, a$variant), ]
  j <- match(lead, r$variants$variant)
  rr <- if (is.na(j)) NULL else fit_variant(r$pheno, r$G[, j], keep_fit = TRUE)
  if (!is.null(rr)) stopifnot(r$variants$alt[j] == drow$alt, r$variants$ref[j] == drow$ref)
  rep_p <- if (is.null(rr)) NA_real_ else rr$row$p
  rep_b <- if (is.null(rr)) NA_real_ else rr$row$beta
  rows[[lead]] <- data.frame(variant = lead, chr = drow$chr, pos = drow$pos, effect_allele = drow$alt,
    discovery_hr = drow$hr, discovery_p = drow$p,
    replication_hr = if (is.null(rr)) NA_real_ else rr$row$hr,
    replication_lower95 = if (is.null(rr)) NA_real_ else rr$row$lower95,
    replication_upper95 = if (is.null(rr)) NA_real_ else rr$row$upper95,
    replication_p = rep_p, replication_threshold = cfg$family_alpha / length(leads),
    replication_status = if (is.null(rr)) 'failed_replication_QC' else rr$row$status,
    same_direction = is.finite(rep_b) && sign(rep_b) == sign(drow$beta),
    replicated = is.finite(rep_p) && rep_p < cfg$family_alpha / length(leads) && sign(rep_b) == sign(drow$beta))
  for (cohort in c('discovery', 'replication')) {
    z <- if (cohort == 'discovery') x else r
    jj <- match(lead, z$variants$variant)
    if (is.na(jj)) next
    ff <- fit_variant(z$pheno, z$G[, jj], keep_fit = TRUE)
    if (!is.null(ff$fit)) {
      ph <- cox.zph(ff$fit)$table
      diagnostics[[paste(lead, cohort)]] <- data.frame(variant = lead, cohort, term = rownames(ph),
                                                     chisq = ph[, 1], df = ph[, 2], p = ph[, 3])
    }
  }
  for (s in unique(x$pheno$site)) {
    ii <- which(x$pheno$site == s); jj <- match(lead, x$variants$variant)
    ss <- fit_variant(x$pheno[ii, ], x$G[ii, jj])$row
    site_rows[[paste(lead, s)]] <- cbind(data.frame(variant = lead, site = s), ss)
  }
}
empty <- data.frame(variant = character(), chr = integer(), pos = integer(), effect_allele = character(),
                     discovery_hr = numeric(), discovery_p = numeric(), replication_hr = numeric(),
                     replication_lower95 = numeric(), replication_upper95 = numeric(), replication_p = numeric(),
                     replication_threshold = numeric(), replication_status = character(), same_direction = logical(), replicated = logical())
write_tsv(if (length(rows)) do.call(rbind, rows) else empty, 'results/replication/lead_results.tsv')
write_tsv(if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(variant = character(), cohort = character(), term = character(), chisq = numeric(), df = numeric(), p = numeric()),
          'results/replication/ph_diagnostics.tsv')
write_tsv(if (length(site_rows)) do.call(rbind, site_rows) else data.frame(variant = character(), site = character()),
          'results/replication/site_sensitivity.tsv')
