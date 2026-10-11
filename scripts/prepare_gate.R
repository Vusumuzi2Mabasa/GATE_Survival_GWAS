# Export QC-passed simulated discovery data without reading simulation truth.
# BED is SNP-major; BIM A1=ALT, A2=REF. Validate every genotype by decoding it.
source('scripts/common.R')
x <- load_analysis('discovery')
out <- 'data/gate'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out, 'discovery')
G <- x$G; v <- x$variants; d <- x$pheno
stopifnot(identical(rownames(G), d$id), identical(colnames(G), v$variant))
write.table(data.frame(d$id, d$id, 0, 0, 0, -9), paste0(prefix, '.fam'),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(data.frame(v$chr, v$variant, 0, v$pos, v$alt, v$ref), paste0(prefix, '.bim'),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
bed <- file(paste0(prefix, '.bed'), 'wb')
writeBin(as.raw(c(108, 27, 1)), bed)
nbyte <- ceiling(nrow(G)/4)
for (j in seq_len(ncol(G))) {
  code <- rep(0L, nbyte * 4L)
  g <- G[, j]
  code[seq_along(g)] <- ifelse(is.na(g), 1L, c(3L, 2L, 0L)[g + 1L])
  packed <- as.integer(crossprod(c(1L, 4L, 16L, 64L), matrix(code, nrow = 4)))
  writeBin(as.raw(packed), bed)
}
close(bed)
bed <- file(paste0(prefix, '.bed'), 'rb')
stopifnot(identical(as.integer(readBin(bed, 'raw', 3)), c(108L,27L,1L)))
for (j in seq_len(ncol(G))) {
  packed <- as.integer(readBin(bed, 'raw', nbyte))
  # vapply columns are bytes; flatten each column in sample order.
  code <- as.vector(vapply(packed, function(b) bitwAnd(bitwShiftR(b, c(0L,2L,4L,6L)), 3L), integer(4)))
  decoded <- c(2, NA, 1, 0)[code[seq_len(nrow(G))] + 1L]
  stopifnot(identical(is.na(decoded), unname(is.na(G[,j]))), all(decoded == G[,j], na.rm = TRUE))
}
close(bed)
# Explicit indicator columns avoid ambiguity in categorical covariate handling.
design <- model.matrix(~ site + ancestry_group, data = d)[, -1, drop = FALSE]
pheno <- cbind(d[, c('id','time','event','null_time','null_event','age10','male','PC1','PC2','PC3')], design)
write_tsv(pheno, paste0(prefix, '_phenotypes.tsv'))
vcf <- file(paste0(prefix, '.vcf'), 'wt')
writeLines(c('##fileformat=VCFv4.2',
  paste0('##contig=<ID=', sort(unique(v$chr)), '>'),
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  paste(c('#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',d$id), collapse='\t')), vcf)
for (j in seq_len(ncol(G))) {
  gt <- c('0/0','0/1','1/1')[G[,j] + 1L]; gt[is.na(gt)] <- './.'
  writeLines(paste(c(v$chr[j],v$pos[j],v$variant[j],v$ref[j],v$alt[j],'.','PASS','.','GT',gt),collapse='\t'),vcf)
}
close(vcf)
writeLines(c(sprintf('Validated %d samples and %d variants by full BED round-trip.',nrow(G),ncol(G)),
 'BED A1=ALT G. VCF uses REF A and ALT G. Missing genotypes preserved.',
 'Event status: 1=event, 0=right censored. Time is years from entry.',
 'VCF is plain text; bgzip and tabix are required before indexed GATE testing.',
 'Sites and ancestry are indicator covariates. This differs from stratified Cox baseline hazards.',
 'Inputs prepared only; this file does not certify a completed GATE analysis.'),
 file.path(out,'validation.txt'))
message('GATE inputs prepared and BED round-trip verified.')
