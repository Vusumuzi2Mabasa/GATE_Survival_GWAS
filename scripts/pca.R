source('scripts/common.R')
cohort <- commandArgs(TRUE)[1]
x <- readRDS(sprintf('results/qc/%s.rds', cohort))
G <- x$G; v <- x$variants
# LD calculation removes between-population frequency differences first.
Z <- ancestry_center(G, x$pheno$ancestry_group)
kept <- integer()
for (j in seq_len(ncol(G))) {
  nearby <- kept[v$chr[kept] == v$chr[j] & abs(v$pos[kept] - v$pos[j]) <= cfg$pca_window_bp]
  r2 <- if (length(nearby)) cor(Z[, j], Z[, nearby, drop = FALSE])^2 else 0
  if (all(is.finite(r2)) && max(r2) < cfg$pca_r2) kept <- c(kept, j)
}
stopifnot(length(kept) > cfg$n_pcs)
# Mean filling is for PCA only; association testing retains missing values.
P <- G[, kept, drop = FALSE]
mu <- colMeans(P, na.rm = TRUE); af <- mu / 2
P <- sweep(P, 2, mu, '-'); P[is.na(P)] <- 0
P <- sweep(P, 2, sqrt(2 * af * (1 - af)), '/')
pc <- prcomp(P, center = FALSE, scale. = FALSE, rank. = cfg$n_pcs)
scores <- scale(pc$x[, 1:cfg$n_pcs, drop = FALSE])
write_tsv(data.frame(id = x$pheno$id, scores), sprintf('results/pca/%s.tsv', cohort))
write_tsv(v[kept, ], sprintf('results/pca/%s_pruned_variants.tsv', cohort))
put_rds(list(pca = pc, means = mu, denominator = sqrt(2 * af * (1 - af)),
             markers = v$variant[kept]), sprintf('results/pca/%s_model.rds', cohort))
png(sprintf('results/pca/%s.png', cohort), 1100, 780, res = 130)
col <- ifelse(x$pheno$ancestry_group == 'Pop_A', '#176C8E', '#C55339')
plot(scores[, 1:2], col = adjustcolor(col, alpha.f = .45), pch = 16, cex = .5,
     xlab = 'PC1 (standardised)', ylab = 'PC2 (standardised)', main = paste('Genotype PCA:', cohort))
legend('topright', c('Simulated Pop A', 'Simulated Pop B'), col = c('#176C8E', '#C55339'), pch = 16, bty = 'n')
dev.off()
