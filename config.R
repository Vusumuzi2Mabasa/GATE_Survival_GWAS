# Prespecified teaching design. No values are estimated from AWI-Gen.
cfg <- list(
  seed = 20261008L, n_discovery = 3600L, n_replication = 2400L,
  chromosomes = 22L, blocks_per_chr = 10L, variants_per_block = 10L,
  sites = 6L, followup_years = 5, weibull_shape = 1.3,
  baseline_scale = 0.025, dropout_rate = 0.035,
  genotype_missing_rate = 0.002,
  sample_missing_max = 0.03, variant_missing_max = 0.02,
  maf_min = 0.05, mac_min = 100L, hwe_min = 1e-6,
  heterozygosity_mad_multiplier = 6,
  pca_r2 = 0.2, pca_window_bp = 250000L, n_pcs = 3L,
  clump_r2 = 0.2, clump_window_bp = 250000L,
  family_alpha = 0.05, conventional_gwas_p = 5e-8,
  causal_chromosomes = c(3L, 10L, 18L),
  causal_block_within_chr = 5L, causal_variant_within_block = 5L,
  causal_hr = c(1.5, 0.67, 1.4)
)
