source('scripts/common.R')
set.seed(cfg$seed)
m <- cfg$chromosomes * cfg$blocks_per_chr * cfg$variants_per_block
nb <- m / cfg$variants_per_block
map <- data.frame(variant = sprintf('sim_rs%05d', seq_len(m)),
                  chr = rep(1:22, each = m / 22),
                  pos = rep(rep(seq_len(cfg$blocks_per_chr) * 200000L,
                                each = cfg$variants_per_block) +
                              rep(seq_len(cfg$variants_per_block) * 1000L, cfg$blocks_per_chr), 22),
                  ref = 'A', alt = 'G')
# Same population frequencies and LD-generating process in independent cohorts.
freq_a <- runif(nb, .15, .40)
freq_b <- pmin(.48, pmax(.07, freq_a + rep(c(-.09, .09), length.out = nb)))
causal <- (cfg$causal_chromosomes - 1L) * (m / 22) +
          (cfg$causal_block_within_chr - 1L) * cfg$variants_per_block + cfg$causal_variant_within_block
truth <- map
truth$block <- rep(seq_len(nb), each = cfg$variants_per_block)
truth$true_beta <- 0; truth$true_beta[causal] <- log(cfg$causal_hr)
truth$causal <- seq_len(m) %in% causal
truth$signal_block <- truth$block %in% truth$block[causal]
fault <- rep('none', m)
fault[1:10] <- 'monomorphic'; fault[11:20] <- 'rare'
fault[21:30] <- 'hwe_excess_heterozygosity'; fault[31:40] <- 'high_missingness'
truth$injected_fault <- fault
write_tsv(map, 'data/raw/variants.tsv')
write_tsv(truth, 'data/truth/variants.tsv')
make_cohort <- function(n, label, seed) {
  set.seed(seed)
  site <- rep(paste0('Site_', 1:6), length.out = n)
  a <- rbinom(n, 1, rep(c(.3, .4, .5, .6, .7, .5), length.out = n))
  G <- matrix(0L, n, m)
  # Diploid haplotypes: two independent ancestral haplotypes per individual.
  # Markers copy a block anchor with probability .93, otherwise redraw.
  for (b in seq_len(nb)) {
    pp <- ifelse(a == 0, freq_a[b], freq_b[b])
    h1 <- rbinom(n, 1, pp); h2 <- rbinom(n, 1, pp)
    for (k in seq_len(cfg$variants_per_block)) {
      j <- (b - 1L) * cfg$variants_per_block + k
      u <- ifelse(runif(n) < .93, h1, rbinom(n, 1, pp))
      v <- ifelse(runif(n) < .93, h2, rbinom(n, 1, pp))
      G[, j] <- u + v
    }
  }
  age <- runif(n, 40, 65); male <- rbinom(n, 1, .45)
  age10 <- (age - 50) / 10
  # Independent causes of baseline hazard: no blood pressure mediator adjustment.
  lp <- log(1.4) * age10 + log(1.15) * male + log(1.35) * a
  baseline <- cfg$baseline_scale * c(.75, .9, 1, 1.1, 1.25, 1.4)[match(site, paste0('Site_', 1:6))]
  genetic_lp <- drop(G[, causal, drop = FALSE] %*% log(cfg$causal_hr))
  T <- (-log(runif(n)) / (baseline * exp(lp + genetic_lp)))^(1 / cfg$weibull_shape)
  T0 <- (-log(runif(n)) / (baseline * exp(lp)))^(1 / cfg$weibull_shape)
  C <- pmin(cfg$followup_years, rexp(n, cfg$dropout_rate))
  C0 <- pmin(cfg$followup_years, rexp(n, cfg$dropout_rate))
  d <- data.frame(id = sprintf('%s_%05d', label, seq_len(n)), site,
                   ancestry_group = ifelse(a == 0, 'Pop_A', 'Pop_B'), age, age10, male,
                   time = pmin(T, C), event = as.integer(T <= C),
                   null_time = pmin(T0, C0), null_event = as.integer(T0 <= C0))
  # Inject outcome-independent errors after generating survival.
  G[, 1:10] <- 0L
  for (j in 11:20) G[, j] <- rbinom(n, 2, .002)
  for (j in 21:30) G[, j] <- rbinom(n, 1, .9)
  G[matrix(runif(n * m) < cfg$genotype_missing_rate, n, m)] <- NA_integer_
  for (j in 31:40) G[runif(n) < .15, j] <- NA_integer_
  for (i in 1:8) G[i, runif(m) < .5] <- NA_integer_
  G[9:12, ] <- 1L
  # Exact duplicate genotype+phenotype records, given different IDs.
  originals <- 101:104
  dd <- d[originals, ]; dd$id <- paste0(dd$id, '_DUP')
  d <- rbind(d, dd); G <- rbind(G, G[originals, , drop = FALSE])
  rownames(G) <- d$id; colnames(G) <- map$variant
  put_rds(list(G = G, pheno = d, variants = map), sprintf('data/raw/%s.rds', label))
  write_tsv(d, sprintf('data/raw/%s_phenotypes.tsv', label))
  write_tsv(data.frame(id = c(d$id[1:12], dd$id),
                       injected_fault = c(rep('high_missingness', 8), rep('high_heterozygosity', 4), rep('duplicate', 4))),
            sprintf('data/truth/%s_samples.tsv', label))
  # Interoperable autosomal VCF. Positions/alleles are fictional, not a reference build.
  con <- gzfile(sprintf('data/raw/%s.vcf.gz', label), 'wt')
  writeLines(c('##fileformat=VCFv4.2', '##source=synthetic_survival_gwas',
               '##reference=SIMULATED_NOT_A_REAL_ASSEMBLY',
               paste0('##contig=<ID=', 1:22, '>'),
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased genotype">',
               paste(c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', d$id), collapse = '\t')), con)
  for (j in seq_len(m)) {
    gt <- c('0/0', '0/1', '1/1')[G[, j] + 1L]; gt[is.na(gt)] <- './.'
    writeLines(paste(c(map$chr[j], map$pos[j], map$variant[j], 'A', 'G', '.', 'PASS', '.', 'GT', gt), collapse = '\t'), con)
  }
  close(con)
  message(label, ': ', nrow(d), ' rows, ', m, ' variants')
}
make_cohort(cfg$n_discovery, 'discovery', cfg$seed + 1L)
make_cohort(cfg$n_replication, 'replication', cfg$seed + 2L)
capture.output(sessionInfo(), file = 'logs/simulation_session.txt')
