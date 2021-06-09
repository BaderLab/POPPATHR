#' Script that recodes PLINK fam file to case/control format by population
#' adapted from SP plink_baseSetup.R (part of GWAS2pathway package)
#'
#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param pop1 (char) character code for the first population (controls)
#' @param pop2 (char) character code for the second population (cases)
#' @param popsF (char) path to file with population information.
#'		Gives the number of samples per population in the dataset.
#' @param phenoPerm (logical) runs GSEA using phenotype permutations
#'		(i.e., randomly re-shuffing case/control labels n times).
#' 		if FALSE, will run GSEA using genotype permutations (see setupGSEArun.R)
#' @param numPerm (integer >=0) number of re-shuffled fam files to write by
#' 		re-shuffling case-control labels in fam file for calcFST.R.
#'    If phenoPerm is FALSE, leave this blank.
#' @param caseCode (integer) value for case samples in fam phenotype
#' 		column
#' @param ctrlCode (integer) value for control samples in fam phenotype
#' 		column
#' @param famName (char) optional - name for fam file. default: case-ctrl.
#'   file extension added.
#' @param setSeed (integer) value for set.seed() before shuffling
#' @return
#' @export

recodeFAM <- function(genoF, pop1, pop2, popsF, recode=TRUE,
											phenoPerm=TRUE, numPerm=100L, setSeed=42L,
											caseCode=2, ctrlCode=1, famName="case-ctrl") {

	if (data == "1KG") { #### NOTE fix use grep here; not using `data` variable anymore ####
		# keep IDs corresponding to population codes
		pops      <- read.delim(popsF, as.is=T, h=T)
		ctrlNames <- pops$Sample.name[pops[,4] == pop1]
		caseNames <- pops$Sample.name[pops[,4] == pop2]
	}
	if (data == "HM3") {
		# keep IDs corresponding to population codes
		pops      <- read.table(popsF, as.is=T, h=F)
		ctrlNames <- pops[,2][pops[,7] == pop1]
		caseNames <- pops[,2][pops[,7] == pop2]
	} else {
		stop("Pipeline only supports HM3 or 1KG SNP genotype data currently.")
	}

	fam_tab <- read.table(sprintf("%s.fam", genoF), h=F)

	if (recode == TRUE) {

		# Write out new fam file w/ case-control coding
		cat("  * Rewriting .fam file with case/control coding...")
		cat(sprintf("In population file:\t %i cases (%s), %i controls (%s)\n",
								length(caseNames), pop2, length(ctrlNames), pop1))

		fam_tab[fam_tab$V2 %in% caseNames, 'V6'] <- caseCode
		fam_tab[fam_tab$V2 %in% ctrlNames, 'V6'] <- ctrlCode

		famDir  <- dirname(genoF)
	  famOut  <- sprintf("%s/%s.fam", famDir, famName)
	  write.table(fam_tab, file=famOut, row=F, quote=F, col=F)

		caseNum <- which(fam_tab$V6 == 2)
		ctrlNum <- which(fam_tab$V6 == 1)

		cat(sprintf(" file written to %s.\n", famOut))
		cat(sprintf("In phenotype file:\t %i cases (%s), %i controls (%s)\n",
								length(caseNum), pop2, length(ctrlNum), pop1))

		# Stop code from continuing if 'pop1/2' individuals from population file
		# are not found in the respective pheno (.fam) file.
		# e.g., GWD (Gambian) individuals are indicated in the 1KG population
		# description file but are not found in the pheno file. This prevents the
		# code from continuing without having any cases / controls.
		if (sum(fam_tab$V6 == 2) == 0)
			 stop(sprintf(paste("No '%s' individuals found in PLINK pheno file.",
			 										"Correct spelling or try another population."), pop2))

		if (sum(fam_tab$V6 == 1) == 0)
	 		 stop(sprintf(paste("No '%s' individuals found in PLINK pheno file.",
			 										"Correct spelling or try another population."), pop1))
	}
	# Now write the permuted files
	if (phenoPerm == TRUE) {

		permDir <- sprintf("%s/perm", outDir)
		if (!file.exists(permDir)) dir.create(permDir)

		fam_orig 	  <- fam_tab
		fam_orig$V6	<- -9

		# redo case/ctrl name vector b/c of errors with 1KG assignment
		caseNames <- as.character(subset(fam_tab$V2, fam_tab$V6 == caseCode))
		ctrlNames <- as.character(subset(fam_tab$V2, fam_tab$V6 == ctrlCode))

		both <- c(caseNames, ctrlNames)
		n <- length(both)

		set.seed(setSeed)
		for (k in 1:numPerm) {
			cat("Shuffling case/control labels in fam file...")

			shuf	  	<- both[sample(n, replace=FALSE)] # shuffle
			shuf_case <- shuf[1:length(caseNames)]
			shuf_ctrl	<- shuf[(length(caseNames)+1):n]
			fam_tab <- fam_orig
			fam_tab[fam_tab$V2%in%shuf_case, 'V6'] <- caseCode
			fam_tab[fam_tab$V2%in%shuf_ctrl, 'V6'] <- ctrlCode

			famOut	<- sprintf("%s/%i.fam", permDir, k)

			write.table(fam_tab, file=famOut, row=F, col=F, quote=F)
			cat(sprintf(" file written to %s/%s.\n", permDir, basename(famOut)))
		}
	} else {
		return(famOut)
	}
}
