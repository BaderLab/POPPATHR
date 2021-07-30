#' Recodes PLINK fam file with population-specific codes
#'
#' @param genotypeFile (char) path to PLINK-formatted SNP genotype data.
#' @param popOne (char) character code for the first population (e.g., CEU).
#' @param popTwo (char) character code for the second population (e.g., YRI).
#' @param populationTable (char) path to file with population information.
#' Gives the number of samples per population in the dataset.
#' @param setSeed (integer) value for set.seed() before shuffling (default=42).
#' @param outputFile (char) name for fam file (default=population-coded).
#' File extension added.
#' @export
recodeFAM <- function(genotypeFile, popOne, popTwo, populationTable,
	setSeed = 42, outputFile = "population-coded") {

	# Read in population table and PLINK fam file
	cat(sprintf("* Reading in population file: %s\n", populationTable))
	cat(sprintf("* Reading in PLINK fam file: %s.fam\n", genotypeFile))
	populations <- read.table(populationTable, as.is=TRUE, header=TRUE)
	famTab <- read.table(sprintf("%s.fam", genotypeFile), header=FALSE)

	# Keep IDs corresponding to specified population codes
	popOneIDs <- populations[,"IID"][populations[,"population"] == popOne]
	popTwoIDs <- populations[,"IID"][populations[,"population"] == popTwo]

	# Write out new fam file w/ case-control coding
	cat(sprintf("* Writing new fam file with population coding (%s and %s genotypes)\n",
		popOne, popTwo))
	famTab[famTab$V2 %in% popOneIDs, "V6"] <- 1
	famTab[famTab$V2 %in% popTwoIDs, "V6"] <- 2

	# Stop code from continuing if popOne/popTwo individuals from population file
	# are not found in the respective phenotype (.fam) file.
	# e.g., GWD (Gambian) individuals are indicated in the 1KG population
	# description file but are not found in the phenotype file.
	if (sum(famTab$V6 == 2) == 0) {
		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 	"Correct spelling or try another population."), popTwo))
	}
	if (sum(famTab$V6 == 1) == 0) {
 		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 	"Correct spelling or try another population."), popOne))
  	}

	# Write out new fam file
	famFolder <- dirname(genotypeFile)
  	famOut <- sprintf("%s/%s.fam", famFolder, outputFile)
  	write.table(famTab, file=famOut, row=FALSE, quote=FALSE, col=FALSE)

	nPopOne <- length(which(famTab$V6 == 1))
	nPopTwo <- length(which(famTab$V6 == 2))
	cat(sprintf("** Found %s %s and %s %s genotypes\n",
		nPopOne, popOne, nPopTwo, popTwo))
	cat(sprintf("* New population-coded fam file written to %s\n", famOut))
}
