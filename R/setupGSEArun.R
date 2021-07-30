#' Sets up and runs GSEA using population-based FST values
#'
#' @param fstFile (char) path to file with SNP-level FST statistics.
#' @param annotationFile (char) path to pathway definitions GMT file.
#' @param snp2geneFile (char) path to file with snp2gene mappings (output of SNP2gene()).
#' @param setPerm (integer) number of permutations to run (default=10000).
#' @param minGene (integer) minimal size of gene set to test (default=10).
#' @param maxGene (integer) maximum size of gene set to test (default=500).
#' @param scoreType (character) defines GSEA score type. Possible options are
#' ("std", "pos", "neg"; default="pos").
#' @param outputFolder (char) path to store output files.
#' @export
setupGSEArun <- function(fstFile, annotationFile, snp2geneFile,
	setPerm=10000, minGene=10, maxGene=500, scoreType=scoreType,
	outputFolder) {

	# Prepare output files
 	resOut <- sprintf("%s/results.txt", outputFolder)
	snpFstOut <- sprintf("%s/snp2gene_fst.txt", outputFolder)
	snpFstTopOut <- sprintf("%s/snp2gene_fstTop.txt", outputFolder)

	# Read in annotations
	anno <- gmtPathways(annotationFile)

	# Prepare data
	cat("* Selecting top mapped SNP (highest FST) per gene\n")
	fst <- read.delim(fstFile, header=TRUE)
	snp <- read.delim(snp2geneFile, header=FALSE)

	#######
	# NOTE testing SNP-gene distance filtering (500kb)
	#snp <- snp[,-which(colnames(snp) %in% "V3")]
	snp <- filter(snp, V3 < 50000)
	snp <- snp[,-which(colnames(snp) %in% "V3")]
	#######

	colnames(snp) <- c("snp", "gene")

	# Combine snp2gene with FST data
	snpFst <- left_join(snp, fst, by="snp")
	snpFstTop <- snpFst[order(snpFst$fst, decreasing=TRUE),]
	snpFstTop <- snpFstTop[!duplicated(snpFstTop$gene),]
	cat(sprintf(" ** %s SNPs remaining from %s total\n", nrow(snpFstTop), nrow(snpFst)))

	# Write out files
	write.table(snpFst, file=snpFstOut, col=TRUE, quote=FALSE, sep="\t")
	write.table(snpFstTop, file=snpFstTopOut, col=TRUE, quote=FALSE, sep="\t")

	# Prepare gene-level data for GSEA
	geneFST <- setNames(snpFstTop$fst, snpFstTop$gene)

  	# Run GSEA
	cat("* Running GSEA\n")
	res <- fgseaSimple(
		pathways=anno,
		stats=geneFST,
		nperm=setPerm,
		minSize=minGene,
		maxSize=maxGene,
		scoreType=scoreType
	)
	res <- res[order(res$NES, decreasing=TRUE),]

	# Write out results file
	cat(sprintf("* Writing out GSEA results file to %s\n", resOut))
	fwrite(res, file=resOut, sep="\t", sep2=c("", " ", ""))
}
