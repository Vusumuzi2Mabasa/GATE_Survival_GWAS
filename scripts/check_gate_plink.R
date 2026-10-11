source('scripts/common.R')
x <- load_analysis('discovery')
f <- read.table('results/gate/plink_input_check.frq',header=TRUE)
f <- f[match(x$variants$variant,f$SNP),]
stopifnot(identical(as.character(f$SNP),x$variants$variant))
af <- ifelse(f$A1 == x$variants$alt, f$MAF, 1-f$MAF)
expected <- colMeans(x$G,na.rm=TRUE)/2
stopifnot(max(abs(af-expected)) < 1e-4,
          all(f$NCHROBS == 2*colSums(!is.na(x$G))))
writeLines(c('PLINK allele frequencies agree with original ALT dosages within output rounding.',
             'Observed allele counts agree for every variant.'),
           'results/gate/plink_validation.txt')
message('Independent PLINK checks passed for all ',nrow(f),' variants.')
