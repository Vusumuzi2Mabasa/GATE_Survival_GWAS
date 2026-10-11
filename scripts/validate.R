source('scripts/common.R')
checks <- data.frame(check = character(), passed = logical(), detail = character())
check <- function(name, condition, detail = '') {
  checks[nrow(checks) + 1L, ] <<- list(name, isTRUE(condition), detail)
}
check('HWE balanced genotypes', hwe_exact(c(rep(0, 25), rep(1, 50), rep(2, 25))) > .9)
check('HWE monomorphic boundary', hwe_exact(rep(0, 100)) == 1)
check('HWE excess heterozygotes', hwe_exact(rep(1, 100)) < 1e-10)
check('HWE allele coding invariant', abs(hwe_exact(c(rep(0, 15), rep(1, 12), rep(2, 3))) -
  hwe_exact(2 - c(rep(0, 15), rep(1, 12), rep(2, 3)))) < 1e-12)
x <- load_analysis('discovery'); r <- load_analysis('replication')
check('Independent cohort IDs', !any(x$pheno$id %in% r$pheno$id))
for (cohort in c('discovery', 'replication')) {
  z <- if (cohort == 'discovery') x else r
  check(paste(cohort, 'sample alignment'), identical(rownames(z$G), z$pheno$id))
  check(paste(cohort, 'variant alignment'), identical(colnames(z$G), z$variants$variant))
  check(paste(cohort, 'positive bounded followup'), all(z$pheno$time > 0 & z$pheno$time <= cfg$followup_years))
  check(paste(cohort, 'finite PCs'), all(is.finite(as.matrix(z$pheno[, paste0('PC', 1:3)]))))
  s <- read_tsv(sprintf('results/qc/%s_samples.tsv', cohort))
  bad <- read_tsv(sprintf('data/truth/%s_samples.tsv', cohort))
  check(paste(cohort, 'injected bad samples excluded'), all(!s$keep[match(bad$id, s$id)]))
  v <- read_tsv(sprintf('results/qc/%s_variants.tsv', cohort))
  truth <- read_tsv('data/truth/variants.tsv')
  check(paste(cohort, 'injected bad variants excluded'), all(!v$keep[truth$injected_fault != 'none']))
}
a <- read_tsv('results/summary/signal_associations.tsv')
a0 <- read_tsv('results/summary/null_associations.tsv')
check('Complete signal scan', setequal(a$variant, x$variants$variant) && !anyDuplicated(a$variant))
check('Complete null scan', setequal(a0$variant, x$variants$variant) && !anyDuplicated(a0$variant))
check('P values in range', all(is.na(a$p) | (a$p >= 0 & a$p <= 1)))
check('Bonferroni includes entire QC family', max(abs(a$bonferroni_p - pmin(1, ifelse(is.na(a$p), 1, a$p) * nrow(a)))) < 1e-12)
check('No failed fits labelled significant', all(!a$scan_significant[a$status != 'ok']))
# Regression check: flipping the counted allele reverses beta, preserving p.
j <- which(a$status == 'ok')[1]
k <- match(a$variant[j], x$variants$variant)
f1 <- fit_variant(x$pheno, x$G[, k])$row
f2 <- fit_variant(x$pheno, 2 - x$G[, k])$row
check('Allele flip changes sign', abs(f1$beta + f2$beta) < 1e-7)
check('Allele flip preserves p', abs(f1$p - f2$p) < 1e-7)
check('Per-variant complete-case N', f1$n == sum(!is.na(x$G[, k])))
# Direct survival::coxph reference, independent of the wrapper formula builder.
dd <- x$pheno; dd$g <- x$G[, k]
reference <- coxph(Surv(time, event) ~ g + age10 + male + PC1 + PC2 + PC3 +
                    strata(site, ancestry_group), data = dd, ties = 'efron')
check('Wrapper agrees with direct Cox fit', abs(unname(coef(reference)['g']) - f1$beta) < 1e-9)
for (file in c('scripts/qc.R', 'scripts/pca.R', 'scripts/scan.R', 'scripts/summarise.R', 'scripts/replicate.R')) {
  check(paste('No ground-truth reads in', file), !any(grepl('data/truth', readLines(file), fixed = TRUE)))
}
rep <- read_tsv('results/replication/lead_results.tsv')
check('Replication restricted to discovery leads', all(rep$variant %in% a$variant[a$scan_significant]))
if (nrow(rep)) check('Replication family correction', max(abs(rep$replication_threshold - cfg$family_alpha / nrow(rep))) < 1e-12)
write_tsv(checks, 'results/validation/checks.tsv')
# Truth is accessed only after association discovery; linked tags are not false positives.
truth <- read_tsv('data/truth/variants.tsv')
joined <- merge(a, truth[, c('variant', 'block', 'causal', 'signal_block', 'true_beta')], by = 'variant')
write_tsv(joined[joined$causal | joined$scan_significant, ], 'results/validation/truth_comparison.tsv')
truth_blocks <- unique(truth$block[truth$causal])
detected_blocks <- unique(joined$block[joined$scan_significant & joined$signal_block])
write_tsv(data.frame(causal_variants = sum(truth$causal), causal_variants_passing_qc = sum(joined$causal),
  causal_variants_significant = sum(joined$causal & joined$scan_significant),
  signal_blocks = length(truth_blocks), detected_signal_blocks = length(detected_blocks),
  significant_variants_outside_signal_blocks = sum(joined$scan_significant & !joined$signal_block),
  null_scan_significant = sum(a0$scan_significant),
  independent_leads = nrow(rep), replicated_leads = sum(rep$replicated)),
  'results/validation/recovery.tsv')
capture.output(sessionInfo(), file = 'logs/analysis_session.txt')
source_files <- c('config.R', 'Snakefile', list.files('scripts', full.names = TRUE, pattern = '\\.[Rr]$'))
write_tsv(data.frame(file = source_files, md5 = unname(tools::md5sum(source_files))), 'results/validation/source_checksums.tsv')
stopifnot(all(checks$passed))
message(nrow(checks), ' computational checks passed. Statistical results are reported, not forced to pass.')
