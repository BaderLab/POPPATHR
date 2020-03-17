#' Recodes PLINK fam file with population-specific codes
#'
#' @param genotype_file (char) path to PLINK-formatted SNP genotype data.
#' @param pop_one (char) character code for the first population (e.g., CEU).
#' @param pop_two (char) character code for the second population (e.g., YRI).
#' @param population_table (char) path to file with population information.
#'		Gives the number of samples per population in the dataset.
#' @param SET_SEED (integer) value for set.seed() before shuffling (default=42).
#' @param out_file (char) name for fam file (default=population-coded).
#'   File extension added.
#'
#' @return none
#' @export
#'

recodeFAM <- function(genotype_file, pop_one, pop_two, population_table,
											SET_SEED=42, out_file="population-coded") {

	# Keep IDs corresponding to population codes
	pops <- read.table(population_table, as.is=TRUE, h=TRUE)
	pop_one_names <- pops[,2][pops[,7] == pop_one]
	pop_two_names <- pops[,2][pops[,7] == pop_two]

	fam_tab <- read.table(sprintf("%s.fam", genotype_file), h=FALSE)

	# Write out new fam file w/ case-control coding
	cat("* Rewriting fam file with population coding...\n")
	cat(sprintf("	In population table: %i cases (%s), %i controls (%s)\n",
		length(pop_two_names), pop_two, length(pop_one_names), pop_one))

	fam_tab[fam_tab$V2 %in% pop_one_names, 'V6'] <- 1
	fam_tab[fam_tab$V2 %in% pop_two_names, 'V6'] <- 2

	# Stop code from continuing if pop_one/pop_two individuals from population file
	# are not found in the respective phenotype (.fam) file.
	# e.g., GWD (Gambian) individuals are indicated in the 1KG population
	# description file but are not found in the phenotype file.
	if (sum(fam_tab$V6 == 2) == 0)
		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 										"Correct spelling or try another population."), pop_two))

	if (sum(fam_tab$V6 == 1) == 0)
 		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 										"Correct spelling or try another population."), pop_one))

	fam_folder <- dirname(genotype_file)
  fam_out <- sprintf("%s/%s.fam", fam_folder, out_file)
  write.table(fam_tab, file=fam_out, row=FALSE, quote=FALSE, col=FALSE)

	caseNum <- which(fam_tab$V6 == 2)
	ctrlNum <- which(fam_tab$V6 == 1)

	cat(sprintf("	In population-coded fam file: %i cases (%s), %i controls (%s)\n",
		length(caseNum), pop_two, length(ctrlNum), pop_one))
	cat(sprintf("* New population-coded file written to %s.\n", fam_out))
}
