source('scripts/common.R')
scenario <- commandArgs(TRUE)[1]
d <- do.call(rbind, lapply(1:22, function(k) read_tsv(sprintf('results/scan/%s/chr%02d.tsv', scenario, k))))
m <- nrow(d)
stopifnot(!anyDuplicated(d$variant))
# Failed fits remain in the prespecified family size. They are not discoveries.
p_for_adjust <- ifelse(is.finite(d$p), d$p, 1)
d$bonferroni_p <- pmin(1, p_for_adjust * m)
d$bh_q <- p.adjust(p_for_adjust, 'BH')
d$scan_significant <- is.finite(d$p) & d$p < cfg$family_alpha / m
d$conventional_genomewide <- is.finite(d$p) & d$p < cfg$conventional_gwas_p
d <- d[order(d$chr, d$pos), ]
write_tsv(d, sprintf('results/summary/%s_associations.tsv', scenario))
ok <- is.finite(d$p)
lambda <- median(qchisq(pmax(d$p[ok], 1e-300), 1, lower.tail = FALSE)) / qchisq(.5, 1)
write_tsv(data.frame(scenario, planned_tests = m, successful_tests = sum(ok),
                     failed_tests = sum(!ok), scan_threshold = cfg$family_alpha / m,
                     significant_variants = sum(d$scan_significant),
                     conventional_hits = sum(d$conventional_genomewide),
                     nominal_p_below_05 = mean(d$p[ok] < .05), lambda_gc = lambda),
          sprintf('results/summary/%s_metrics.tsv', scenario))
png(sprintf('results/summary/%s_manhattan.png', scenario), 1500, 700, res = 140)
par(mar = c(5, 5, 4, 2))
xx <- seq_len(m); yy <- -log10(safe_p(d$p))
plot(xx, yy, pch = 16, cex = .55, col = c('#176C8E', '#71AAB5')[1 + d$chr %% 2], xaxt = 'n',
     xlab = 'Simulated chromosome', ylab = expression(-log[10](p)),
     ylim = c(0, max(yy, -log10(cfg$conventional_gwas_p), na.rm = TRUE) * 1.12),
     main = paste('Survival association scan:', scenario))
axis(1, at = tapply(xx, d$chr, mean), labels = 1:22)
abline(h = -log10(cfg$family_alpha / m), col = '#C55339', lty = 2)
abline(h = -log10(cfg$conventional_gwas_p), col = 'grey35', lty = 3)
legend('topright', c('Bonferroni for simulated scan', 'Conventional 5e-8 reference'),
       lty = c(2, 3), col = c('#C55339', 'grey35'), bty = 'n', cex = .8)
dev.off()
png(sprintf('results/summary/%s_qq.png', scenario), 850, 800, res = 120)
observed <- -log10(safe_p(sort(d$p[ok]))); expected <- -log10(ppoints(sum(ok)))
plot(expected, observed, pch = 16, cex = .55, col = '#176C8E',
     xlab = 'Expected -log10(p)', ylab = 'Observed -log10(p)', main = paste('QQ plot:', scenario))
abline(0, 1, col = '#C55339', lty = 2)
mtext(sprintf('Descriptive lambda GC = %.3f; LD makes tests dependent', lambda), cex = .8)
dev.off()
