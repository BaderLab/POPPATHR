#!/usr/bin/Rscript

#' writes subset of gmt with feature-selected pathways
#'
#' @param groupList (list) names are input nets, values are member entities
#' e.g. pathways and genes
#' @param cutoff (integer) only features with scores >= cutoff will be
#' included
#' @param units (char) units that were present in the data. Not all units
#' in groupList may be in the data. This function ensures that the units
#' written out are limited to those for which data were available
#' e.g. for a pathway, units are the genes that were present in the
#' dataset. Genes that were not assayed would not be included in the output.
#' To include all genes (not recommended unless you know what you are doing)
#' , set to "*".
#' @param outFile (char) where output should be written.
#' @export
#-------------------------------------------------------------------------------
# Input files
dataDir <- "/media/catherine/DATAPART1/Data/PatientNetworks/Catherine"
outDir  <- sprintf("%s/GWAS_GSEA/PNC", dataDir)

snp2geneF  <- sprintf("%s/snp2gene/snp2gene.txt", outDir)
pathList <- sprintf("%s/anno/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
		    dataDir)
outFile	<- sprintf("%s/updated_anno.gmt", outDir)
#-------------------------------------------------------------------------------

writeTopGroups <- function(pathList, snp2geneFile, outFile) {

	require(PatientClassifier)
	pathways <- readPathways(pathList, MIN_SIZE=0, MAX_SIZE=20000)
	cat("* reading snp file\n")
	snp2gene <- read.delim(snp2geneF, h=F, as.is=T, sep="\t")
	#genes <- subset(snp2gene, select=V2)
	genes <- snp2gene[,2]

	nGenes <- length(unique(unlist(pathways)))
	cat(sprintf("%i genes are mapped to SNPs", length(genes)))
	cat(sprintf(" %i pathways contain %i unique genes\n", length(pathways), nGenes))

	system(sprintf("cat /dev/null > %s",outFile))
	t0 <- Sys.time()

	num_genes <- numeric()
		for (p in names(pathways)) {
				g <- intersect(pathways[[p]],genes)
				num_genes <- c(num_genes, length(g))
				cat(sprintf("%s: %i genes intersect\n", p, length(g)))
				cat(sprintf("%s\t%s\t%s\n",p,p,paste(g,collapse="\t")),
										file=outFile,append=TRUE)

			warning("don't do this!")
			stop("you need to now get the top genes, and take thm out of the gmt file")

	}

	t1 <- Sys.time()
	print(t1-t0)
	print(summary(num_genes))

	cat(sprintf("GMT written to %s\n", outFile))
}

writeTopGroups(pathList, statFile, outFile)
