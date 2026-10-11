source('scripts/common.R')
library(ggplot2)
library(gtsummary)
library(flextable)
library(htmltools)
dir.create('report/plots', recursive = TRUE, showWarnings = FALSE)
dir.create('report/tables', recursive = TRUE, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(),
  plot.title = element_text(face = 'bold', color = '#123B4A'),
  plot.caption = element_text(hjust = 0, color = '#52646D'), legend.position = 'bottom'))
save_plot <- function(p, name, width = 11, height = 5.4) {
  ggsave(file.path('report/plots', name), p, width = width, height = height, dpi = 160, bg = 'white')
}
x <- load_analysis('discovery'); r <- load_analysis('replication')
sig <- read_tsv('results/summary/signal_associations.tsv')
null <- read_tsv('results/summary/null_associations.tsv')
metrics <- rbind(read_tsv('results/summary/signal_metrics.tsv'), read_tsv('results/summary/null_metrics.tsv'))
qc <- rbind(read_tsv('results/qc/discovery_summary.tsv'), read_tsv('results/qc/replication_summary.tsv'))
rep <- read_tsv('results/replication/lead_results.tsv')
recovery <- read_tsv('results/validation/recovery.tsv')
checks <- read_tsv('results/validation/checks.tsv')
ph <- read_tsv('results/replication/ph_diagnostics.tsv')
calibration <- read_tsv('results/validation/null_calibration_summary.tsv')
null_reps <- read_tsv('results/validation/null_repetitions.tsv')
d <- rbind(transform(x$pheno, cohort = 'Discovery'), transform(r$pheno, cohort = 'Replication'))
d$sex <- factor(d$male, 0:1, c('Female', 'Male'))
d$event_status <- factor(d$event, 0:1, c('Censored', 'Event'))
cohort_table <- tbl_summary(d, by = cohort,
  include = c(age, sex, site, ancestry_group, time, event_status),
  label = list(age ~ 'Age (years)', sex ~ 'Simulated sex', site ~ 'Study site',
               ancestry_group ~ 'Simulated ancestry group', time ~ 'Observed follow-up (years)',
               event_status ~ 'Signal-outcome status'),
  statistic = list(all_continuous() ~ '{median} ({p25}, {p75})', all_categorical() ~ '{n} ({p}%)'),
  missing = 'ifany')
cohort_ft <- as_flex_table(cohort_table)
style_ft <- function(ft) {
  ft <- theme_booktabs(ft)
  ft <- bg(ft, part = 'header', bg = '#123B4A')
  ft <- color(ft, part = 'header', color = 'white')
  ft <- fontsize(ft, size = 10, part = 'all')
  ft <- font(ft, fontname = 'Arial', part = 'all')
  autofit(ft)
}
make_ft <- function(d) {
  if (!nrow(d)) d <- data.frame(Result = 'No rows for this analysis.')
  style_ft(flextable(d))
}
cohort_ft <- style_ft(cohort_ft)
qc_display <- data.frame(Measure = c('Initial sample records', 'Samples retained', 'Initial markers',
  'Markers retained', 'Signal-outcome events', 'Null-outcome events'),
  Discovery = as.numeric(unlist(qc[qc$cohort == 'discovery', -1])),
  Replication = as.numeric(unlist(qc[qc$cohort == 'replication', -1])))
qc_ft <- make_ft(qc_display)
metric_labels <- c('Planned variant tests', 'Successful fits', 'Failed fits', 'Scan significance threshold',
  'Scan-significant variants', 'Variants below 5e-8', 'Fraction with p < 0.05', 'Lambda GC')
format_metrics <- function(s) {
  z <- metrics[metrics$scenario == s, ]
  c(format(unlist(z[c('planned_tests', 'successful_tests', 'failed_tests')]), trim = TRUE),
    format(z$scan_threshold, scientific = TRUE, digits = 3),
    z$significant_variants, z$conventional_hits, sprintf('%.3f', z$nominal_p_below_05), sprintf('%.3f', z$lambda_gc))
}
metric_display <- data.frame(Measure = metric_labels, Signal = format_metrics('signal'), Null = format_metrics('null'))
metrics_ft <- make_ft(metric_display)
rep_display <- if (nrow(rep)) data.frame(
  Variant = rep$variant, Chromosome = rep$chr,
  'Discovery HR' = round(rep$discovery_hr, 3),
  'Discovery p' = format.pval(rep$discovery_p, digits = 3, eps = 1e-300),
  'Replication HR (95% CI)' = sprintf('%.3f (%.3f to %.3f)', rep$replication_hr,
    rep$replication_lower95, rep$replication_upper95),
  'Replication p' = format.pval(rep$replication_p, digits = 3, eps = 1e-300),
  Replicated = ifelse(rep$replicated, 'Yes', 'No'), check.names = FALSE) else data.frame()
rep_ft <- make_ft(rep_display)
truth_ft <- make_ft(data.frame(Measure = c('Planted direct-effect variants', 'Direct-effect variants passing QC',
  'Direct-effect variants significant', 'Planted signal blocks', 'Detected signal blocks',
  'Significant variants outside signal blocks', 'Significant variants in full null scan',
  'Discovery lead regions', 'Replicated lead variants'), Value = as.numeric(unlist(recovery[1, ]))))
checks_ft <- make_ft(data.frame(Check = checks$check, Result = ifelse(checks$passed, 'Pass', 'Fail')))
saveRDS(cohort_table, 'report/tables/cohort_gtsummary.rds')
saveRDS(list(cohort = cohort_ft, qc = qc_ft, association = metrics_ft,
             replication = rep_ft, recovery = truth_ft, checks = checks_ft), 'report/tables/flextables.rds')
# Direct R HTML export avoids a Pandoc dependency on this computer.
tables_page <- tags$html(lang = 'en', tags$head(tags$title('Synthetic survival GWAS tables'),
  tags$meta(charset = 'utf-8'), tags$meta(name = 'viewport', content = 'width=device-width, initial-scale=1'),
  tags$style(HTML('body{font-family:Arial,sans-serif;padding:24px;color:#243d48}section{overflow-x:auto;margin-bottom:32px}h1,h2{color:#123b4a}'))),
  tags$body(h1('Synthetic survival GWAS tables'),
    tags$section(h2('Cohort characteristics'), htmltools_value(cohort_ft)),
    tags$section(h2('QC flow'), htmltools_value(qc_ft)),
    tags$section(h2('Association scans'), htmltools_value(metrics_ft)),
    tags$section(h2('Independent replication'), htmltools_value(rep_ft))))
save_html(tables_page, file = 'report/tables/tables.html', libdir = 'lib')
for (scenario in c('signal', 'null')) {
  a <- if (scenario == 'signal') sig else null
  a$x <- seq_len(nrow(a)); a$y <- -log10(safe_p(a$p))
  ticks <- aggregate(x ~ chr, a, mean)
  plot <- ggplot(a, aes(x, y, color = factor(chr %% 2))) +
    geom_point(size = .9, alpha = .8, na.rm = TRUE) +
    scale_color_manual(values = c('#176C8E', '#7CB4C0'), guide = 'none') +
    scale_x_continuous(breaks = ticks$x, labels = ticks$chr) +
    geom_hline(yintercept = -log10(cfg$family_alpha / nrow(a)), linetype = 2, color = '#C55339') +
    geom_hline(yintercept = -log10(cfg$conventional_gwas_p), linetype = 3, color = '#52646D') +
    labs(title = paste('Survival GWAS training scan:', scenario), x = 'Simulated chromosome', y = expression(-log[10](p)),
      subtitle = 'Common variants; adjusted Cox proportional hazards models',
      caption = sprintf('Dashed: simulated-scan Bonferroni %.3g. Dotted: 5e-8 reference. All data are synthetic.', cfg$family_alpha / nrow(a)))
  save_plot(plot, paste0(scenario, '_manhattan.png'))
}
qq <- do.call(rbind, lapply(c('signal', 'null'), function(s) {
  a <- if (s == 'signal') sig else null
  p <- sort(a$p[is.finite(a$p)])
  data.frame(scenario = s, expected = -log10(ppoints(length(p))), observed = -log10(safe_p(p)))
}))
save_plot(ggplot(qq, aes(expected, observed)) + geom_abline(slope = 1, intercept = 0, color = '#C55339', linetype = 2) +
  geom_point(size = 1, color = '#176C8E') + facet_wrap(~scenario) +
  labs(title = 'Association calibration and signal', x = 'Expected -log10(p)', y = 'Observed -log10(p)',
       caption = 'The separate null outcome has no genetic effects. One scan is not a formal type-I-error study.'), 'qq.png')
save_plot(ggplot(null_reps, aes(replicate, nominal_rate)) +
  geom_hline(yintercept = .05, linetype = 2, color = '#C55339') +
  geom_line(color = '#176C8E') + geom_point(color = '#176C8E', size = 2) +
  scale_x_continuous(breaks = seq(1, nrow(null_reps), by = 2)) +
  labs(title = 'Repeated null outcomes on a fixed LD-pruned marker panel', x = 'Independent outcome repetition',
       y = 'Fraction of marker p-values below 0.05',
       caption = 'Dashed line: nominal 5%. Repetitions regenerate outcomes, not genotypes. Dependence within each panel is retained.'), 'null_calibration.png')
save_plot(ggplot(d, aes(PC1, PC2, color = ancestry_group)) + geom_point(size = .6, alpha = .5) +
  facet_wrap(~cohort) + scale_color_manual(values = c('#176C8E', '#C55339')) +
  labs(title = 'Population structure after genotype QC', color = 'Simulated ancestry',
       caption = 'PCA uses ancestry-centred LD pruning, followed by globally standardised genotypes; cohorts are fitted separately.'), 'pca.png')
km <- survfit(Surv(time, event) ~ cohort, data = d)
ks <- summary(km, censored = TRUE)
km_data <- data.frame(time = ks$time, survival = ks$surv, cohort = sub('cohort=', '', as.character(ks$strata)))
km_data <- rbind(data.frame(time = 0, survival = 1, cohort = c('Discovery', 'Replication')), km_data)
km_data <- km_data[order(km_data$cohort, km_data$time), ]
save_plot(ggplot(km_data, aes(time, survival, color = cohort)) + geom_step(linewidth = .9) +
  scale_color_manual(values = c('#176C8E', '#C55339')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = 'Event-free survival in the synthetic cohorts', x = 'Years since cohort entry',
       y = 'Kaplan-Meier survival probability', color = 'Cohort',
       caption = 'Descriptive cohort curves; no genotype selection or causal comparison is implied.'), 'survival.png')
kr <- summary(km, times = c(0, 1, 3, 5))
risk_table <- data.frame(Cohort = sub('cohort=', '', as.character(kr$strata)), Year = kr$time, At_risk = kr$n.risk)
write_tsv(risk_table, 'report/tables/numbers_at_risk.tsv')
risk_ft <- make_ft(risk_table)
if (nrow(rep)) {
  save_plot(ggplot(rep, aes(replication_hr, reorder(variant, replication_hr))) +
    geom_vline(xintercept = 1, linetype = 2, color = '#52646D') +
    geom_segment(aes(x = replication_lower95, xend = replication_upper95, yend = reorder(variant, replication_hr)), color = '#176C8E') +
    geom_point(aes(shape = replicated), size = 3, color = '#176C8E') + scale_x_log10() +
    labs(title = 'Independent replication of discovery lead variants', x = 'Hazard ratio per effect allele (log scale)', y = NULL,
         caption = 'Intervals are nominal 95% CIs; replication decisions use multiplicity-adjusted p-value criteria.'), 'replication.png', height = 4.5)
}
ft_block <- function(ft) div(class = 'table-scroll', htmltools_value(ft))
figure <- function(file, alt) tags$figure(tags$img(src = paste0('plots/', file), alt = alt))
sources <- list(
  c('Cox models in the survival package', 'https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html'),
  c('gwasurvivr: genome-wide survival analysis', 'https://bioconductor.org/packages/release/bioc/html/gwasurvivr.html'),
  c('GATE: scalable survival GWAS with relatedness', 'https://www.nature.com/articles/s41467-022-32885-x'),
  c('PLINK: ancestry-aware HWE quality control', 'https://www.cog-genomics.org/plink/2.0/tutorials/qc1'),
  c('PLINK: genotype principal components', 'https://www.cog-genomics.org/plink/2.0/strat'),
  c('gtsummary to flextable', 'https://www.danieldsjoberg.com/gtsummary/reference/as_flex_table.html'))
page <- tags$html(lang = 'en', tags$head(tags$title('Simulated survival GWAS | R analysis report'),
  tags$meta(charset = 'utf-8'), tags$meta(name = 'viewport', content = 'width=device-width, initial-scale=1'),
  tags$style(HTML('body{margin:0;background:#f1f5f5;color:#243d48;font-family:Arial,sans-serif;line-height:1.65}main{max-width:1180px;margin:32px auto;background:white;padding:40px;border-radius:12px}h1,h2{color:#123b4a;line-height:1.25}h1{font-size:36px}h2{margin-top:42px;border-top:1px solid #dce6e7;padding-top:24px}.label{color:#176c8e;font-weight:bold;letter-spacing:1px}.note{background:#edf5f7;border-left:4px solid #176c8e;padding:16px 22px}.table-scroll{overflow-x:auto;margin:20px 0}figure{margin:24px 0}img{max-width:100%;height:auto}a{color:#126681}code{background:#edf3f3;padding:2px 5px}li{margin:8px 0}@media(max-width:700px){main{margin:0;padding:20px}h1{font-size:28px}}'))),
  tags$body(tags$main(
    div(class = 'label', 'SYNTHETIC DATA / LOCAL TRAINING PIPELINE'),
    h1('Survival GWAS from simulation to replication'),
    p('A reproducible R analysis coordinated by Snakemake. Tables use gtsummary and flextable; figures use ggplot2. No AWI-Gen or other participant data are used.'),
    div(class = 'note', p(sprintf('%s participants passed discovery QC; %s observed events. %s variants were scanned across 22 simulated chromosomes.',
      format(nrow(x$pheno), big.mark = ','), sum(x$pheno$event), nrow(sig))),
      p(sprintf('%d discovery lead regions; %d independently replicated lead variants. %d of %d planted signal blocks detected.',
        nrow(rep), sum(rep$replicated), recovery$detected_signal_blocks, recovery$signal_blocks))),
    h2('1. Prespecified design and estimand'),
    p('The outcome is time from disease-free cohort entry to a simulated diabetes-like event, with exact event times, independent dropout and administrative censoring at five years. Baseline age and sex affect the hazard. Three variants have conditional hazard ratios of 1.50, 0.67 and 1.40; these deliberately large effects make the teaching example informative.'),
    p('Two fictional ancestry groups occur within six sites. Both ancestry and site affect baseline hazards. Independent haplotypes produce correlated variants in LD blocks. Unrelated people are generated; copied records and technical faults are injected afterward. The groups are not models of named African populations.'),
    p('Each scan fits one variant dosage plus age, sex and three genotype PCs, with separate baseline hazards for each site-by-ancestry stratum. Effects are associations, not causal claims. Single-variant estimates need not equal conditional generating effects when other causal variants are omitted.'),
    h2('2. Cohorts and quality control'), ft_block(cohort_ft),
    figure('survival.png', 'Kaplan-Meier event-free survival curves for the synthetic cohorts'), ft_block(risk_ft), ft_block(qc_ft),
    p('Sample QC removes excessive missingness, within-ancestry heterozygosity outliers, exact copied records and invalid outcomes. Variant QC checks missingness, minor-allele frequency, minor-allele count and exact Hardy-Weinberg equilibrium within ancestry groups. Filters are fixed in config.R and do not use outcome associations.'),
    p('Sex-chromosome checks, imputation, reference harmonisation and general kinship estimation are not represented by this autosomal hard-call simulation. Exact duplicate detection is not a substitute for KING or another validated relatedness method on real data.'),
    figure('pca.png', 'PCA scatter plots for independent discovery and replication cohorts'),
    h2('3. Genome-wide association scans'), ft_block(metrics_ft),
    p('The primary test is a two-sided common-variant Cox Wald test. Missing genotypes are excluded separately for each variant. Failed or warned fits are flagged and retained in the multiple-testing family size. Efron handling is used for tied event times.'),
    p('The primary training-scan threshold is 0.05 divided by all QC-passing variants. The conventional 5e-8 threshold is displayed separately. A reduced synthetic scan is not equivalent to a densely imputed human GWAS. BH q-values are supplied as exploratory outputs, not substituted for the primary threshold.'),
    figure('signal_manhattan.png', 'Manhattan plot for the simulated genetic signal outcome'),
    figure('null_manhattan.png', 'Manhattan plot for the independent null outcome'),
    figure('qq.png', 'QQ plots of signal and null association p values'),
    p('The null outcome retains age, sex, site and ancestry effects but has no genetic effect within strata. Its false positives and lambda GC are descriptive checks. LD makes marker tests dependent; one simulated cohort cannot establish power or repeated-sampling error control. No genomic-control rescaling is applied.'),
    figure('null_calibration.png', 'Nominal rejection fractions across independent null-outcome repetitions'),
    p(sprintf('Across %d independent null-outcome repetitions on %d LD-pruned markers, the mean nominal rejection fraction was %.3f (Monte Carlo SE %.3f; approximate 95%% interval %.3f to %.3f). Mean lambda GC was %.3f.',
      calibration$repetitions, calibration$markers_per_repetition, calibration$mean_nominal_rate,
      calibration$monte_carlo_se, calibration$nominal_rate_lower95, calibration$nominal_rate_upper95, calibration$mean_lambda_gc)),
    p(sprintf('The smaller pruned panel had family-wise false-positive frequency %.3f (exact binomial 95%% interval %.3f to %.3f), using its own 0.05/m threshold. This is not an estimate of full-GWAS family-wise error. These repetitions provide only a coarse check, conditional on this fixed genotype sample.',
      calibration$pruned_panel_fwer, calibration$pruned_panel_fwer_lower95, calibration$pruned_panel_fwer_upper95)),
    h2('4. Independent replication and diagnostics'), ft_block(rep_ft),
    p('Discovery leads are selected using a prespecified distance and within-ancestry LD threshold. All discovery leads count toward the replication correction, including leads unavailable after replication QC. Replication requires a two-sided p below 0.05 divided by the number of leads and matching effect direction. Alleles are explicitly checked.'),
    if (nrow(rep)) figure('replication.png', 'Replication hazard ratios and confidence intervals') else NULL,
    p(sprintf('Schoenfeld-residual diagnostics were recorded for %d cohort-by-variant global tests; %d had nominal p < 0.05. These selected-model diagnostics are exploratory, not a guarantee of proportional hazards.',
      sum(ph$term == 'GLOBAL'), sum(ph$p[ph$term == 'GLOBAL'] < .05, na.rm = TRUE))),
    p('Site-specific sensitivity estimates are saved in results/replication/site_sensitivity.tsv. Sparse sites may fail to estimate an effect; failures are reported. Replication is an independent draw from the same simulation model, not external validation in a new real-world population.'),
    h2('5. Recovery of known truth and computational checks'), ft_block(truth_ft), ft_block(checks_ft),
    p('Known causal labels are used only here and in validation, after discovery and replication. Correlated tags in a planted signal block are not counted as false discoveries merely because their direct generating coefficient is zero. LD clumping is not fine-mapping.'),
    h2('6. Reproducibility and next experiments'),
    p('Run ./run_pipeline.ps1 from PowerShell, or bash run_snakemake.sh from WSL in this project. All source, configuration, logs, results and the isolated Snakemake environment live on the Toshiba drive. Source checksums and R session information are retained.'),
    tags$ul(tags$li('Repeat the simulation across independent seeds to quantify power, family-wise error and effect-estimate bias.'),
      tags$li('Add interval censoring, competing risks or relatedness only with a matching analysis method.'),
      tags$li('Use scalable validated tools such as GATE/SPACox for appropriate large-cohort settings; do not simply extrapolate this full-fit R scanner to millions of variants.'),
      tags$li('Keep prediction modelling (for example LASSO-Cox with nested validation) separate from association discovery.')),
    h2('Methods and software sources'), tags$ul(lapply(sources, function(s) tags$li(tags$a(href = s[2], s[1])))),
    p(tags$a(href = '../results/summary/signal_associations.tsv', 'Full discovery results'), ' | ',
      tags$a(href = '../results/replication/lead_results.tsv', 'Replication results'), ' | ',
      tags$a(href = 'tables/tables.html', 'Open the formatted tables separately'), ' | ',
      tags$a(href = '../README.md', 'Project documentation'))
  )))
save_html(page, file = 'report/index.html', libdir = 'lib')
capture.output(sessionInfo(), file = 'logs/report_session.txt')
write_tsv(data.frame(package = c('gtsummary', 'flextable', 'ggplot2', 'htmltools', 'survival'),
  version = vapply(c('gtsummary', 'flextable', 'ggplot2', 'htmltools', 'survival'), function(p) as.character(packageVersion(p)), character(1))),
  'report/tables/package_versions.tsv')
message('R report generated: report/index.html')
