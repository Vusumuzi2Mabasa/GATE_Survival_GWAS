library(GATE)
args <- commandArgs(trailingOnly=TRUE)
outcome <- if(length(args)) args[1] else 'signal'
stopifnot(outcome %in% c('signal','null'))
prefix <- paste0('results/gate/discovery_',outcome)
e <- new.env()
load(paste0(prefix,'.rda'), envir=e)
stopifnot(exists('modglmm',envir=e,inherits=FALSE))
m <- e$modglmm
d <- read.delim('data/gate/discovery_phenotypes.tsv')
vr <- scan(paste0(prefix,'.varianceRatio.txt'),quiet=TRUE)
stopifnot(isTRUE(m$converged),identical(m$traitType,'survival'),
          length(m$sampleID)==nrow(d),setequal(as.character(m$sampleID),d$id),
          all(is.finite(m$theta)),all(m$theta>=0),
          length(vr)>0,all(is.finite(vr)),all(vr>0))
write.table(data.frame(samples=nrow(d),events=sum(d[[if(outcome=='signal') 'event' else 'null_event']]),converged=m$converged,
                       variance_ratio=vr[1],theta=paste(m$theta,collapse=',')),
            paste0(prefix,'_model_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
message('Null-model convergence, participant IDs and variance ratio passed checks.')
