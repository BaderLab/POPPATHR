#' Maps SNPs to their nearest genes
#'
#' NOTE: updated GenomicRanges to version 1.30.1 on 01/30/2018
#' Requires ignore.strand=TRUE param to properly run distanceToNearest()
#' given unknown strand assignment for SNP array-based genotyping info
#'
#' @param in_file (char) path to PLINK bim file (only CHR, SNP, and BP
#' 		columns considered).
#' @param refgene_file (char) path to refseq table with header.
#' @param marg (integer) region upstream and downstream of [txStart,txEnd] (default=0).
#' @param out_file (char) path to write snp-to-gene mapping file.
#'
#' @return none
#' @export
#'

SNP2gene <- function(in_file, refgene_file, marg=0, out_file) {

	# Read in SNP table
	cat("* Reading SNP table\n")
	snps <- fread(in_file, h=FALSE, data.table=FALSE)
	snp_GR  <- GRanges(paste("chr", snps[,1], sep=""),
					  				 IRanges(snps[,4], snps[,4]),
					  			 	 name=snps[,2])

	refGene <- read.delim(refgene_file, sep="\t", h=TRUE, as.is=TRUE)
	gene_GR	<- GRanges(refGene[,"chrom"],
										 IRanges(refGene[,"txStart"] + 1 - marg,
										 refGene[,"txEnd"] + marg),
										 strand=refGene[,"strand"],
										 name=refGene[,"name2"])

	start_GR	<- resize(gene_GR, fix="start", width=1L)
	end_GR		<- resize(gene_GR, fix="end", width=1L)

	cat("* Computing distances\n")
	# SNPs inside the gene domain
	d0		<- distanceToNearest(snp_GR, gene_GR, ignore.strand=TRUE)
	dbit	<- d0@elementMetadata$distance
	d_in	<- data.frame(queryHits=d0@from,
											subjectHits=d0@to,
											distance=dbit)

	idx		<- which(dbit == 0)
	out1	<- cbind(snp_GR$name[d_in$queryHits[idx]],
					 			 gene_GR$name[d_in$subjectHits[idx]],
					 		 	 d_in$distance[idx])

	# SNPs outside gene domain
	idx		<- which(dbit > 0) # not in a gene
	snp2_GR <- snp_GR[idx]
	cat(sprintf("%i SNPs not inside gene domain\n", length(idx)))

	rm(snp_GR, snps, idx)

	cat("* Computing distance to domain starts\n")
	d1	 <- distanceToNearest(snp2_GR, start_GR, ignore.strand=TRUE)
	dbit <- d1@elementMetadata$distance
	d_start <- cbind(d1@from, d1@to, dbit)
	colnames(d_start)[1:2] <- c("queryHits", "subjectHits")

	cat("* Computing distance to domain ends\n")
	d2	 <- distanceToNearest(snp2_GR, end_GR, ignore.strand=TRUE)
	dbit <- d2@elementMetadata$distance
	d_end	<- cbind(d2@from, d2@to, dbit)
	colnames(d_end)[1:2] <- c("queryHits","subjectHits")

	if (all.equal(d_start[,1], d_end[,1]) != TRUE) {
		cat("d_start and d_end have different indexing. You need some other ")
		cat("way to include all snps\n")
	}

	d_both	<- merge(d_start, d_end, by.x="queryHits", by.y="queryHits")

	# Start is closer than end
	idx		<- which(d_both$dbit.x < d_both$dbit.y);
	out2	<- cbind(snp2_GR$name[d_both$queryHits[idx]],
					 			 start_GR$name[d_both$subjectHits.x[idx]],
					 		 	 d_both$dbit.x[idx])
	idx		<- setdiff(1:nrow(d_both), idx)
	out3	<- cbind(snp2_GR$name[d_both$queryHits[idx]],
					 			 end_GR$name[d_both$subjectHits.y[idx]],
					 		 	 d_both$dbit.y[idx])
	out <- rbind(out1, out2, out3)

	cat(sprintf("* Writing output to file %s\n", out_file))
	options(scipen=10)  # prevent conversion of numbers to sci notation
	write.table(out, file=out_file, sep="\t", col=FALSE, row=FALSE, quote=FALSE)
}
