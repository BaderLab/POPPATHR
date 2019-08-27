#' Recodes PLINK fam file to case/control format by population
#' adapted from SP plink_baseSetup.R (part of GWAS2pathway package)
#'
#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param pop1 (char) character code for the first population (controls).
#' @param pop2 (char) character code for the second population (cases).
#' @param popsF (char) path to file with population information.
#'		Gives the number of samples per population in the dataset.
#' @param caseCode (integer) value for case samples in fam phenotype column.
#' @param ctrlCode (integer) value for control samples in fam phenotype column.
#' @param outF (char) optional - name for fam file (default=case-ctrl).
#'   file extension added.
#' @param setSeed (integer) value for set.seed() before shuffling (default=42).
#'
#' @return none
#' @export
#'

recodeFAM <- function(genoF, pop1, pop2, popsF,
											setSeed=42L, caseCode=2, ctrlCode=1,
											outF="case-ctrl") {
	# Keep IDs corresponding to population codes
	pops      <- read.table(popsF, as.is=TRUE, h=TRUE)
	ctrlNames <- pops[,2][pops[,7] == pop1]
	caseNames <- pops[,2][pops[,7] == pop2]

	fam_tab <- read.table(sprintf("%s.fam", genoF), h=FALSE)

	# Write out new fam file w/ case-control coding
	cat("*Rewriting .fam file with case/control coding...\n")
	cat(sprintf("  In population file: %i cases (%s), %i controls (%s)\n",
		length(caseNames), pop2, length(ctrlNames), pop1))

	fam_tab[fam_tab$V2 %in% caseNames, 'V6'] <- caseCode
	fam_tab[fam_tab$V2 %in% ctrlNames, 'V6'] <- ctrlCode

	# Stop code from continuing if 'pop1/2' individuals from population file
	# are not found in the respective pheno (.fam) file.
	# e.g., GWD (Gambian) individuals are indicated in the 1KG population
	# description file but are not found in the pheno file. This prevents the
	# code from continuing without having any cases / controls.
	if (sum(fam_tab$V6 == 2) == 0)
		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 										"Correct spelling or try another population."), pop2))

	if (sum(fam_tab$V6 == 1) == 0)
 		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 										"Correct spelling or try another population."), pop1))

	famDir <- dirname(genoF)
  famOut <- sprintf("%s/%s.fam", famDir, outF)
  write.table(fam_tab, file=famOut, row=F, quote=F, col=F)

	caseNum <- which(fam_tab$V6 == 2)
	ctrlNum <- which(fam_tab$V6 == 1)

	cat(sprintf("  In phenotype fam file: %i cases (%s), %i controls (%s)\n",
		length(caseNum), pop2, length(ctrlNum), pop1))
	cat(sprintf("*File written to %s.\n", famOut))
}
