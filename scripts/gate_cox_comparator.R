# Match the pilot's fixed covariates, common baseline and mean-imputed dosages.
source('scripts/common.R')
x <- load_analysis('discovery')
dir.create('results/gate', recursive = TRUE, showWarnings = FALSE)
idx <- which(x$variants$chr == 3)
result <- lapply(idx, function(j) {
  d <- x$pheno
  d$g <- x$G[,j]
  d$g[is.na(d$g)] <- mean(d$g, na.rm = TRUE)
  fit <- coxph(Surv(time,event) ~ g + age10 + male + PC1 + PC2 + PC3 + site + ancestry_group,
               data=d, ties='efron')
  b <- coef(fit)['g']; se <- sqrt(vcov(fit)['g','g'])
  data.frame(variant=x$variants$variant[j],chr=3,pos=x$variants$pos[j],
             effect_allele=x$variants$alt[j], beta=unname(b), se=se,
             hr=exp(unname(b)), p=2*pnorm(-abs(b/se)), n=nrow(d),events=sum(d$event))
})
write_tsv(do.call(rbind,result),'results/gate/chr3_matched_cox.tsv')
message('Matched Cox comparator completed for ',length(idx),' variants.')
