#' Maps SNPs to their nearest genes and filters SNPs that lie
#' >20,000 bp away from its nearest mapped gene (for use in setupGSEArun.R)
#' (adapted from map_SNP2gene from GWAS2pathway package)
#'
#' @param snpFile (char) path to PLINK .bim file (only CHR, SNP, and BP
#' 		columns considered
#' @param geneFile (char) path to refseq table with header
#' @param marg (integer) region upstream and downstream of [txStart,txEnd]
#' @param filterSNP (logical) if TRUE filters SNP list based on distance
#'		maximum distance from nearest gene
#' @param filterDist (integer) distance cutoff to filter SNP-gene pair
#' 		(i.e., don't include SNP if it lies more than 20KB away from gene)
#' @param outFile (char) path to write snp2gene mapping
#' @export
#' @return No value. Side effect of writing file

require(data.table)
require(GenomicRanges)

args     <- commandArgs(TRUE)
snpFile  <- args[1]
geneFile <- args[2]
outFile  <- args[3]

marg       <- 0L
filterSNP  <- FALSE
filterDist <- 20000L

cat("* Reading SNP table\n")
snps <- fread(snpFile, h=F, data.table=F) #fread much quicker than read.table
snp_GR	<- GRanges(paste("chr", snps[,1], sep=""),
									 IRanges(snps[,4], snps[,4]),
									 name=snps[,2])

refGene <- read.delim(geneFile, sep="\t", h=T, as.is=T)
gene_GR	<- GRanges(refGene[,"chrom"],
									 IRanges(refGene[,"txStart"] + 1 - marg,
									 refGene[,"txEnd"] + marg),
									 strand=refGene[,"strand"],
									 name=refGene[,"name2"])

start_GR	<- resize(gene_GR, fix="start", width=1L)
end_GR		<- resize(gene_GR, fix="end", width=1L)

cat("* Computing distances\n")
# snps inside the gene domain
d0		<- distanceToNearest(snp_GR, gene_GR)
dbit	<- d0@elementMetadata@listData$distance
d_in	<- data.frame(
			 queryHits=d0@queryHits,
			 subjectHits=d0@subjectHits,
			 distance=dbit)

idx		<- which(dbit == 0)
out1	<- cbind(snp_GR$name[d_in$queryHits[idx]],
							 gene_GR$name[d_in$subjectHits[idx]],
							 d_in$distance[idx])

# snps outside gene domain
idx		<- which(dbit > 0) # not in a gene
snp2_GR <- snp_GR[idx]
cat(sprintf("%i SNPs not inside gene domain\n", length(idx)))

rm(snp_GR, snps, idx)

cat("* Computing distance to domain starts\n")
d1	 <- distanceToNearest(snp2_GR, start_GR)
dbit <- d1@elementMetadata@listData$distance
d_start <- cbind(d1@queryHits, d1@subjectHits, dbit)
colnames(d_start)[1:2] <- c("queryHits","subjectHits")

cat("* Computing distance to domain ends\n")
d2	 <- distanceToNearest(snp2_GR, end_GR)
dbit <- d2@elementMetadata@listData$distance
d_end	<- cbind(d2@queryHits, d2@subjectHits, dbit)
colnames(d_end)[1:2] <- c("queryHits","subjectHits")

if (all.equal(d_start[,1], d_end[,1]) != TRUE) {
	cat("d_start and d_end have different indexing. You need some other ")
	cat("way to include all snps\n")
}

d_both	<- merge(d_start, d_end, by.x="queryHits", by.y="queryHits")

# start is closer than end
idx		<- which(d_both$dbit.x < d_both$dbit.y);
out2	<- cbind(snp2_GR$name[d_both$queryHits[idx]],
							 start_GR$name[d_both$subjectHits.x[idx]],
							 d_both$dbit.x[idx]
				)
idx		<- setdiff(1:nrow(d_both), idx)
out3	<- cbind(snp2_GR$name[d_both$queryHits[idx]],
							 end_GR$name[d_both$subjectHits.y[idx]],
							 d_both$dbit.y[idx]
				)

if (filterSNP == TRUE) {
	# Filter out SNPs that lie far away from nearest gene.
	# Only relevant to 'out2' and 'out3' since 'out1' includes SNPs directly
	# within gene domains (i.e., distance b/w SNP and gene = 0).
	out2_new <- out2[as.numeric(out2[,3]) <= filterDist, ]
	out3_new <- out3[as.numeric(out3[,3]) <= filterDist, ]
	out <- rbind(out1, out2_new, out3_new)
	cat(sprintf(paste("Filtered out %i SNPs from %i total SNPs that are >%ibp",
										"away from its nearest mapped gene (%i new total)\n"),
				sum( (nrow(out2)-nrow(out2_new)) + (nrow(out3)-nrow(out3_new)) ),
				sum( nrow(out1) + nrow(out2) + nrow(out3) ),
				filterDist, nrow(out)))

	cat("* Writing to output file\n")
	options(scipen=10) #don't convert numbers to sci notation
	write.table(out, file=outFile, sep="\t", col=F, row=F, quote=F)
} else {
	out <- rbind(out1, out2, out3)
	cat("* Writing to output file\n")
	options(scipen=10) #don't convert numbers to sci notation
	write.table(out, file=outFile, sep="\t", col=F, row=F, quote=F)
}
