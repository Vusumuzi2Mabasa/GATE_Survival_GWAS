# Repeated null-outcome checks conditional on one fixed, QC-passing genotype set.
# This is a calibration experiment, not a genome-wide discovery or power study.
source('scripts/common.R')
x <- load_analysis('discovery')
markers <- read_tsv('results/pca/discovery_pruned_variants.tsv')$variant
jj <- match(markers, colnames(x$G)); stopifnot(!anyNA(jj))
args <- commandArgs(TRUE)
B <- if (length(args)) as.integer(args[1]) else 20L
stopifnot(is.finite(B), B >= 2L)
results <- vector('list', B)
for (b in seq_len(B)) {
  set.seed(cfg$seed + 1000L + b)
  d <- x$pheno; n <- nrow(d)
  a <- as.integer(d$ancestry_group == 'Pop_B')
  lp <- log(1.4) * d$age10 + log(1.15) * d$male + log(1.35) * a
  baseline <- cfg$baseline_scale * c(.75, .9, 1, 1.1, 1.25, 1.4)[match(d$site, paste0('Site_', 1:6))]
  T <- (-log(runif(n)) / (baseline * exp(lp)))^(1 / cfg$weibull_shape)
  C <- pmin(cfg$followup_years, rexp(n, cfg$dropout_rate))
  d$null_time <- pmin(T, C); d$null_event <- as.integer(T <= C)
  p <- vapply(jj, function(j) fit_variant(d, x$G[, j], 'null')$row$p, numeric(1))
  ok <- is.finite(p)
  results[[b]] <- data.frame(replicate = b, seed = cfg$seed + 1000L + b,
    events = sum(d$null_event), planned_tests = length(jj), successful_tests = sum(ok),
    nominal_rejections = sum(p[ok] < .05), nominal_rate = sum(p[ok] < .05) / length(jj),
    panel_bonferroni_threshold = cfg$family_alpha / length(jj),
    panel_any_false_positive = any(p[ok] < cfg$family_alpha / length(jj)),
    lambda_gc = median(qchisq(pmax(p[ok], 1e-300), 1, lower.tail = FALSE)) / qchisq(.5, 1))
  message('Null repetition ', b, '/', B, ': nominal rate ', round(results[[b]]$nominal_rate, 3))
}
z <- do.call(rbind, results)
mean_rate <- mean(z$nominal_rate)
mcse <- sd(z$nominal_rate) / sqrt(B)
ci <- pmax(0, pmin(1, mean_rate + c(-1, 1) * qt(.975, B - 1) * mcse))
fw <- binom.test(sum(z$panel_any_false_positive), B)$conf.int
write_tsv(z, 'results/validation/null_repetitions.tsv')
write_tsv(data.frame(repetitions = B, markers_per_repetition = length(jj),
  mean_nominal_rate = mean_rate, monte_carlo_se = mcse,
  nominal_rate_lower95 = ci[1], nominal_rate_upper95 = ci[2],
  mean_lambda_gc = mean(z$lambda_gc),
  pruned_panel_fwer = mean(z$panel_any_false_positive),
  pruned_panel_fwer_lower95 = fw[1], pruned_panel_fwer_upper95 = fw[2],
  failed_tests = sum(z$planned_tests - z$successful_tests)),
  'results/validation/null_calibration_summary.tsv')
