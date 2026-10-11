# R-authored supplement to the original Cox report.
library(officer)
library(flextable)
library(ggplot2)
read <- function(p) read.delim(p,check.names=FALSE)
s <- read('results/gate/summary.tsv')
a <- read('results/gate/signal_all.tsv')
n <- read('results/gate/null_all.tsv')
c <- read('results/gate/chr3_comparison.tsv')
m <- read('results/gate/discovery_signal_model_checks.tsv')
mn <- read('results/gate/discovery_null_model_checks.tsv')
dir.create('report/gate',recursive=TRUE,showWarnings=FALSE)
a$x <- seq_len(nrow(a))
centres <- aggregate(x~CHR,a,mean)
p <- ggplot(a,aes(x,-log10(p.value),colour=factor(CHR%%2)))+geom_point(size=.8)+
  geom_hline(yintercept=-log10(s$threshold[1]),linetype=2)+
  scale_x_continuous(breaks=centres$x,labels=centres$CHR)+
  scale_colour_manual(values=c('#17365D','#5B9BD5'),guide='none')+
  labs(x='Chromosome',y=expression(-log[10](p)),title='GATE survival association scan')+
  theme_minimal(base_size=12)
ggsave('report/gate/manhattan.png',p,width=7,height=3.7,dpi=200)
qq <- do.call(rbind,lapply(list(signal=a$p.value,null=n$p.value),function(v)
  data.frame(expected=-log10(ppoints(length(v))),observed=-log10(sort(v)))))
qq$outcome <- rep(c('Signal outcome','Null outcome'),each=nrow(a))
p <- ggplot(qq,aes(expected,observed,colour=outcome))+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_point(size=.9)+scale_colour_manual(values=c('#666666','#17365D'))+
  labs(x='Expected negative log10 p',y='Observed negative log10 p',colour=NULL)+
  theme_minimal(base_size=12)+theme(legend.position='bottom')
ggsave('report/gate/qq.png',p,width=7,height=3.7,dpi=200)
p <- ggplot(c,aes(cox_beta,gate_beta))+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_point(colour='#17365D')+labs(x='Matched Cox log hazard ratio',y='GATE effect estimate')+
  theme_minimal(base_size=12)
ggsave('report/gate/comparison.png',p,width=7,height=3.7,dpi=200)
doc <- read_docx()
doc <- body_set_default_section(doc,prop_section(page_size=page_size(8.5,11),
  page_margins=page_mar(top=.65,bottom=.65,left=.75,right=.75)))
for(style in c('Normal','Title','heading 1')) {
  id <- styles_info(doc)$style_id[styles_info(doc)$style_name==style][1]
  doc <- docx_set_paragraph_style(doc,style_id=id,style_name=style,
    fp_t=fp_text(font.family='Arial',font.size=if(style=='Title') 22 else if(style=='heading 1') 15 else 11,
      color='#000000',bold=style!='Normal'),
    fp_p=fp_par(padding.bottom=7,keep_with_next=style!='Normal'))
}
para <- function(t,style='Normal') doc <<- body_add_par(doc,t,style)
heading <- function(t) para(t,'heading 1')
image <- function(p) doc <<- body_add_img(doc,p,width=7,height=3.7)
table <- function(d,widths) {
  ft <- flextable(d)
  ft <- font(ft,fontname='Arial',part='all'); ft <- fontsize(ft,size=10,part='all')
  ft <- padding(ft,padding=5,part='all'); ft <- bg(ft,bg='#17365D',part='header')
  ft <- color(ft,color='white',part='header')
  ft <- bg(ft,bg=rep(c('white','#F0F4F8'),length.out=nrow(d)),part='body')
  ft <- border_outer(ft,fp_border(color='#D9D9D9',width=.5),part='all')
  ft <- border_inner(ft,fp_border(color='#D9D9D9',width=.5),part='all')
  ft <- width(ft,width=widths); ft <- valign(ft,valign='center',part='all')
  ft <- set_table_properties(ft,layout='fixed',opts_word=list(split=FALSE,repeat_headers=TRUE))
  doc <<- body_add_flextable(doc,ft); para('')
}
para('GATE simulated survival GWAS results','Title')
para('R generated analysis supplement for Vusumuzi Mabasa')
para('This report documents the completed synthetic survival GWAS workflow, from data generation and quality control to GATE installation, model fitting, association testing and validation. Detailed methodology and the execution record precede the results. The analysis uses fictional data and does not use AWI Gen participants or measurements.')
para('Methods expanded on 11 October 2026. The analysis described below was completed before this documentation update; no association models were rerun for the update.')
source('scripts/gate_methodology.R',local=TRUE)
doc <- body_add_break(doc)
heading('Results')
para(sprintf('GATE 0.45 completed scans of %s variants for both a signal outcome and a null outcome in %s simulated participants. The signal scan identified %d markers at the panel Bonferroni threshold; the null scan identified %d.',
  format(s$variants[1],big.mark=','),format(m$samples,big.mark=','),s$significant[1],s$significant[2]))
heading('Validation and scan results')
para(sprintf('Both survival null models converged. The signal outcome had %d events and the null outcome had %d. All %d output checks passed for each scan, including participant counts, variant identities, positions, effect alleles, allele frequencies, finite estimates and required SPA convergence.',m$events,mn$events,s$checks_passed[1]))
table(data.frame(Outcome=s$outcome,Variants=s$variants,Significant=s$significant,
  `Lambda GC`=round(s$lambda,3),`SPA tests`=s$spa_required,check.names=FALSE),c(1.4,1.2,1.4,1.4,1.6))
para(sprintf('The threshold is %.3g, calculated as 0.05 divided by %d variants. The conventional 5e-8 threshold yielded %d signal markers and %d null markers. Lambda GC is descriptive here; one null dataset does not establish calibration.',s$threshold[1],s$variants[1],s$significant_5e8[1],s$significant_5e8[2]))
para('Both lambda values are below one, indicating deflation. This requires further investigation with repeated null simulations and a LOCO comparison before production use.')
heading('Planted variants')
z <- a[match(c('sim_rs00245','sim_rs00945','sim_rs01745'),a$SNPID),]
table(data.frame(Variant=z$SNPID,`Effect allele`=z$Allele2,`GATE HR`=round(exp(z$BETA),3),
  `P value`=format(z$p.value,digits=3,scientific=TRUE),check.names=FALSE),c(2.3,1.3,1.5,1.9))
para('Variants within a region are correlated through simulated linkage disequilibrium. Significant markers therefore do not represent independent discoveries. Effect estimates are exponentiated GATE coefficients, with G as the effect allele.')
doc <- body_add_break(doc)
heading('Association patterns')
image('report/gate/manhattan.png')
para('The dashed line is the panel Bonferroni threshold. Results use a genetic relationship matrix fitted without leaving out the tested chromosome.')
image('report/gate/qq.png')
para('The signal outcome departs from the null line at its strongest associations. Interpret the null plot together with its descriptive lambda and the limitations below.')
doc <- body_add_break(doc)
heading('Comparison with Cox and interpretation')
image('report/gate/comparison.png')
lead <- c[c$variant=='sim_rs00245',]
para(sprintf('For the planted chromosome 3 variant sim_rs00245, GATE estimated HR %.3f with p %.3g; the matched Cox model estimated HR %.3f with p %.3g. The 100 variant comparison uses identical covariates and mean imputation. GATE uses a frailty mixed model and SPA score testing; Cox uses a fixed effects model and Wald testing, so their estimates and p values need not coincide.',lead$gate_hr,lead$gate_p,lead$cox_hr,lead$cox_p))
heading('Methods and limits')
para('Events are coded 1 and right censoring 0. Observed times are in years and are not binned. Covariates are age per ten years, sex, three principal components, and explicit site and ancestry indicators. Genotypes passed the original pipeline QC and PLINK export checks. CSI indexing is required for the regional VCF reader used here.')
para(sprintf('The estimated variance ratios were %.4f for the signal model and %.4f for the null model. Variance components were %s and %s. GATE used 30 markers to estimate each ratio.',m$variance_ratio,mn$variance_ratio,m$theta,mn$theta))
para('This is a non LOCO implementation trial with a small marker panel and unrelated simulated people. Including tested regions in the relationship matrix can absorb their signal. A production analysis needs suitable relationship markers, a validated LOCO strategy, relatedness and ancestry assessment, repeated null simulations, and independent replication. This run does not demonstrate biobank scale performance or calibrated extreme tail error rates. GATE is not a Bayesian posterior sampling model.')
para('Sources: github.com/weizhou0/GATE and the GATE paper, Nature Communications 2022, doi 10.1038/s41467-022-32885-x. Source commit 962f1028120c73e24ec34f29e02eeae3bdc3bfb3. Build compatibility changes and package versions are recorded in logs/gate.')
print(doc,target='report/GATE_Survival_GWAS_Results.docx')
writeLines(capture.output(sessionInfo()),'report/gate/R_session.txt')
