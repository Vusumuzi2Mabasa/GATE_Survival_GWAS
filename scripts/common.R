options(stringsAsFactors = FALSE, warn = 1)
library(survival)
source('config.R')
stopifnot(cfg$chromosomes == 22L, cfg$sites == 6L, cfg$n_pcs == 3L)
dir.create('results', showWarnings = FALSE)
dir.create('logs', showWarnings = FALSE)
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
}
read_tsv <- function(path) read.delim(path, check.names = FALSE)
put_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path, compress = FALSE)
}
# Wigginton-style conditional exact HWE test, two-sided probability ordering.
# Applied within the two simulated homogeneous ancestry groups, never pooled.
hwe_exact <- function(g) {
  g <- g[!is.na(g)]
  obs_hom1 <- sum(g == 0); obs_hets <- sum(g == 1); obs_hom2 <- sum(g == 2)
  rare <- 2 * min(obs_hom1, obs_hom2) + obs_hets
  n <- length(g)
  if (n == 0 || rare == 0) return(1)
  probs <- numeric(rare + 1L)
  mid <- floor(rare * (2 * n - rare) / (2 * n))
  if ((mid %% 2) != (rare %% 2)) mid <- mid + 1L
  probs[mid + 1L] <- 1
  h <- mid; r <- (rare - mid) / 2; c <- n - h - r
  while (h >= 2) {
    probs[h - 1L] <- probs[h + 1L] * h * (h - 1) / (4 * (r + 1) * (c + 1))
    h <- h - 2; r <- r + 1; c <- c + 1
  }
  h <- mid; r <- (rare - mid) / 2; c <- n - h - r
  while (h <= rare - 2) {
    probs[h + 3L] <- probs[h + 1L] * 4 * r * c / ((h + 2) * (h + 1))
    h <- h + 2; r <- r - 1; c <- c - 1
  }
  probs <- probs / sum(probs)
  min(1, sum(probs[probs <= probs[obs_hets + 1L] + 1e-12]))
}
ancestry_center <- function(G, groups) {
  Z <- matrix(0, nrow(G), ncol(G), dimnames = dimnames(G))
  for (a in unique(groups)) {
    ii <- which(groups == a)
    v <- G[ii, , drop = FALSE]
    mu <- colMeans(v, na.rm = TRUE); mu[!is.finite(mu)] <- 0
    v <- sweep(v, 2, mu, '-')
    v[is.na(v)] <- 0
    Z[ii, ] <- v
  }
  Z
}
model_formula <- function(with_g = TRUE, scenario = 'signal') {
  response <- if (scenario == 'null') 'Surv(null_time, null_event)' else 'Surv(time, event)'
  as.formula(paste(response, '~', if (with_g) 'g +' else '',
                   'age10 + male + PC1 + PC2 + PC3 + strata(site, ancestry_group)'))
}
# Common-variant Cox scan; failed or warned fits receive no inferential p-value.
# Missing genotypes use per-variant complete cases, never mean-filled association data.
fit_variant <- function(d, g, scenario = 'signal', keep_fit = FALSE) {
  take <- !is.na(g)
  z <- d[take, , drop = FALSE]; z$g <- g[take]
  ev <- if (scenario == 'null') z$null_event else z$event
  af <- mean(z$g) / 2
  row <- data.frame(n = nrow(z), events = sum(ev), effect_allele_frequency = af,
                    mac = 2 * nrow(z) * min(af, 1 - af),
                    event_carriers = sum(z$g[ev == 1] > 0),
                    beta = NA_real_, se = NA_real_, hr = NA_real_,
                    lower95 = NA_real_, upper95 = NA_real_, p = NA_real_,
                    status = 'not_fitted', warning = '')
  if (sum(ev) < 30 || length(unique(z$g)) < 2) {
    row$status <- 'insufficient_information'; return(list(row = row, fit = NULL))
  }
  warns <- character()
  fit <- tryCatch(withCallingHandlers(
    coxph(model_formula(TRUE, scenario), data = z, ties = 'efron',
          singular.ok = FALSE, control = coxph.control(iter.max = 40),
          x = keep_fit, model = keep_fit),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart('muffleWarning') }),
    error = function(e) e)
  if (inherits(fit, 'error')) {
    row$status <- 'error'; row$warning <- conditionMessage(fit)
    return(list(row = row, fit = NULL))
  }
  if (length(warns)) {
    row$status <- 'warning'; row$warning <- paste(unique(warns), collapse = ' | ')
    return(list(row = row, fit = NULL))
  }
  b <- coef(fit)['g']; se <- sqrt(vcov(fit)['g', 'g'])
  if (!is.finite(b) || !is.finite(se) || se <= 0) {
    row$status <- 'nonfinite'; return(list(row = row, fit = NULL))
  }
  row$beta <- b; row$se <- se; row$hr <- exp(b)
  row$lower95 <- exp(b - 1.96 * se); row$upper95 <- exp(b + 1.96 * se)
  row$p <- 2 * pnorm(abs(b / se), lower.tail = FALSE); row$status <- 'ok'
  list(row = row, fit = if (keep_fit) fit else NULL)
}
load_analysis <- function(cohort) {
  x <- readRDS(sprintf('results/qc/%s.rds', cohort))
  pcs <- read_tsv(sprintf('results/pca/%s.tsv', cohort))
  stopifnot(identical(x$pheno$id, pcs$id), identical(rownames(x$G), pcs$id))
  x$pheno <- cbind(x$pheno, pcs[, paste0('PC', 1:cfg$n_pcs), drop = FALSE])
  x
}
safe_p <- function(p) pmax(p, .Machine$double.xmin)
escape_html <- function(s) {
  s <- gsub('&', '&amp;', as.character(s), fixed = TRUE)
  s <- gsub('<', '&lt;', s, fixed = TRUE)
  gsub('>', '&gt;', s, fixed = TRUE)
}
html_table <- function(d) {
  if (!nrow(d)) return('<p>No rows.</p>')
  for (j in seq_along(d)) if (is.numeric(d[[j]])) d[[j]] <- format(signif(d[[j]], 4), trim = TRUE)
  paste0('<table><thead><tr>', paste0('<th>', escape_html(names(d)), '</th>', collapse = ''),
         '</tr></thead><tbody>', paste(apply(d, 1, function(r)
           paste0('<tr>', paste0('<td>', escape_html(r), '</td>', collapse = ''), '</tr>')), collapse = ''),
         '</tbody></table>')
}
