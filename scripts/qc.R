source('scripts/common.R')
cohort <- commandArgs(TRUE)[1]
x <- readRDS(sprintf('data/raw/%s.rds', cohort))
G <- x$G; d <- x$pheno; v <- x$variants
stopifnot(identical(rownames(G), d$id), identical(colnames(G), v$variant),
          !anyDuplicated(d$id), !anyDuplicated(v$variant),
          all(na.omit(as.vector(G)) %in% 0:2))
miss <- rowMeans(is.na(G)); het <- rowMeans(G == 1, na.rm = TRUE)
het_outlier <- rep(FALSE, nrow(G))
for (a in unique(d$ancestry_group)) {
  ii <- which(d$ancestry_group == a & miss <= cfg$sample_missing_max)
  center <- median(het[ii]); spread <- mad(het[ii])
  het_outlier[ii] <- abs(het[ii] - center) > cfg$heterozygosity_mad_multiplier * max(spread, .005)
}
# Exact copied records only. This is not a general kinship estimator.
fingerprint <- apply(G, 1, paste, collapse = ',')
dup <- duplicated(fingerprint)
valid_pheno <- complete.cases(d) & d$time > 0 & d$time <= cfg$followup_years &
               d$event %in% 0:1 & d$null_time > 0 & d$null_time <= cfg$followup_years & d$null_event %in% 0:1
keep <- miss <= cfg$sample_missing_max & !het_outlier & !dup & valid_pheno
sample_qc <- data.frame(id = d$id, missing_rate = miss, heterozygosity = het,
                        heterozygosity_outlier = het_outlier, exact_duplicate = dup,
                        valid_phenotype = valid_pheno, keep = keep)
write_tsv(sample_qc, sprintf('results/qc/%s_samples.tsv', cohort))
G <- G[keep, , drop = FALSE]; d <- d[keep, , drop = FALSE]
miss_v <- colMeans(is.na(G)); af <- colMeans(G, na.rm = TRUE) / 2
maf <- pmin(af, 1 - af); mac <- 2 * colSums(!is.na(G)) * maf
hp <- sapply(unique(d$ancestry_group), function(a)
             apply(G[d$ancestry_group == a, , drop = FALSE], 2, hwe_exact))
min_hwe <- apply(hp, 1, min)
keep_v <- is.finite(maf) & miss_v <= cfg$variant_missing_max & maf >= cfg$maf_min &
          mac >= cfg$mac_min & min_hwe >= cfg$hwe_min
variant_qc <- cbind(v, data.frame(missing_rate = miss_v, effect_allele_frequency = af,
                                 maf = maf, mac = mac, min_within_ancestry_hwe_p = min_hwe, keep = keep_v))
write_tsv(variant_qc, sprintf('results/qc/%s_variants.tsv', cohort))
put_rds(list(G = G[, keep_v, drop = FALSE], pheno = d, variants = v[keep_v, ]),
        sprintf('results/qc/%s.rds', cohort))
write_tsv(data.frame(cohort, initial_samples = nrow(x$G), retained_samples = nrow(d),
                     initial_variants = ncol(x$G), retained_variants = sum(keep_v),
                     events = sum(d$event), null_events = sum(d$null_event)),
          sprintf('results/qc/%s_summary.tsv', cohort))
