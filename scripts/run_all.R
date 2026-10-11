# Dependency-free sequential alternative to Snakemake, using the same R stages.
# Run from project root: Rscript --vanilla scripts/run_all.R
exe <- file.path(R.home('bin'), 'Rscript')
if (.Platform$OS.type == 'windows') exe <- paste0(exe, '.exe')
dir.create('logs', showWarnings = FALSE)
run <- function(script, args = character()) {
  log <- file.path('logs', paste0(basename(script), '_', paste(args, collapse = '_'), '.log'))
  status <- system2(exe, c('--vanilla', shQuote(script), vapply(args, shQuote, character(1))), stdout = log, stderr = log)
  if (status != 0) stop('Stage failed: ', script, '; see ', log)
}
run('scripts/simulate.R')
for (cohort in c('discovery', 'replication')) {
  run('scripts/qc.R', cohort); run('scripts/pca.R', cohort)
}
for (scenario in c('signal', 'null')) {
  for (chr in 1:22) run('scripts/scan.R', c(scenario, sprintf('%02d', chr)))
  run('scripts/summarise.R', scenario)
}
run('scripts/replicate.R'); run('scripts/null_calibration.R', '20')
run('scripts/validate.R'); run('scripts/report.R')
run('scripts/report_word.R')
