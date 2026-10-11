source('scripts/common.R')
x <- load_analysis('discovery')
summaries <- list()
for(outcome in c('signal','null')) {
  a <- do.call(rbind,lapply(1:22,function(chr) read.table(
    sprintf('results/gate/chr%d_%s_full.tsv',chr,outcome),header=TRUE,check.names=FALSE)))
  j <- match(a$SNPID,x$variants$variant)
  expected_af <- colMeans(x$G,na.rm=TRUE)/2
  spa_needed <- abs(a$Tstat/sqrt(a$varT)) > 2
  checks <- c(variant_count=nrow(a)==nrow(x$variants), unique_ids=!anyDuplicated(a$SNPID),
    ids=all(!is.na(j)), participants=all(a$N==nrow(x$pheno)),
    chromosome=all(a$CHR==x$variants$chr[j]),position=all(a$POS==x$variants$pos[j]),
    reference=all(a$Allele1==x$variants$ref[j]),effect_allele=all(a$Allele2==x$variants$alt[j]),
    allele_frequency=all(abs(a$AF_Allele2-expected_af[j])<1e-6),
    finite_effects=all(is.finite(a$BETA)&is.finite(a$SE)&a$SE>0),
    valid_p=all(is.finite(a$p.value)&a$p.value>0&a$p.value<=1),
    positive_variance=all(is.finite(a$varT)&a$varT>0),
    spa_convergence=all(a$Is.SPA.converge[spa_needed]==1))
  write_tsv(data.frame(check=names(checks),passed=unname(checks)),paste0('results/gate/',outcome,'_checks.tsv'))
  stopifnot(all(checks))
  a <- a[order(a$CHR,a$POS),]
  write_tsv(a,paste0('results/gate/',outcome,'_all.tsv'))
  summaries[[outcome]] <- data.frame(outcome=outcome,variants=nrow(a),
    threshold=.05/nrow(a),significant=sum(a$p.value<.05/nrow(a)),
    significant_5e8=sum(a$p.value<5e-8),lambda=median(qchisq(a$p.value,1,lower.tail=FALSE))/qchisq(.5,1),
    spa_required=sum(spa_needed),spa_converged=sum(a$Is.SPA.converge==1),checks_passed=sum(checks))
}
write_tsv(do.call(rbind,summaries),'results/gate/summary.tsv')
a <- read_tsv('results/gate/signal_all.tsv')
c <- read_tsv('results/gate/chr3_matched_cox.tsv')
b <- a[a$CHR==3,]; j <- match(b$SNPID,c$variant)
stopifnot(!anyNA(j), all(b$Allele2==c$effect_allele[j]))
comp <- data.frame(variant=b$SNPID,gate_beta=b$BETA,gate_se=b$SE,gate_hr=exp(b$BETA),
  gate_p=b$p.value,cox_beta=c$beta[j],cox_hr=c$hr[j],cox_p=c$p[j])
write_tsv(comp,'results/gate/chr3_comparison.tsv')
print(do.call(rbind,summaries))
print(comp[comp$variant=='sim_rs00245',])
