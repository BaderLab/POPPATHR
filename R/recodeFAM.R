#' Recodes PLINK fam file with population-specific codes
#'
#' @param genotype_file (char) path to PLINK-formatted SNP genotype data.
#' @param pop_one (char) character code for the first population (e.g., CEU).
#' @param pop_two (char) character code for the second population (e.g., YRI).
#' @param population_table (char) path to file with population information.
#'		Gives the number of samples per population in the dataset.
#' @param SET_SEED (integer) value for set.seed() before shuffling (default=42).
#' @param output_file (char) name for fam file (default=population-coded).
#'   File extension added.
#'
#' @return none
#' @export
#'

recodeFAM <- function(genotype_file, pop_one, pop_two, population_table,
											SET_SEED=42, output_file="population-coded") {

	# Read in population table and PLINK fam file
	cat(sprintf("* Reading in population file: %s\n", population_table))
	cat(sprintf("* Reading in PLINK fam file: %s.fam\n", genotype_file))
	populations <- read.table(population_table, as.is=TRUE, h=TRUE)
	fam_tab <- read.table(sprintf("%s.fam", genotype_file), h=FALSE)

	# Keep IDs corresponding to specified population codes
	pop_one_ids <- populations[,"IID"][populations[,"population"] == pop_one]
	pop_two_ids <- populations[,"IID"][populations[,"population"] == pop_two]

	# Write out new fam file w/ case-control coding
	cat(sprintf("* Writing new fam file with population coding (%s and %s genotypes)\n",
		pop_one, pop_two))

	fam_tab[fam_tab$V2 %in% pop_one_ids, "V6"] <- 1
	fam_tab[fam_tab$V2 %in% pop_two_ids, "V6"] <- 2

	# Stop code from continuing if pop_one/pop_two individuals from population file
	# are not found in the respective phenotype (.fam) file.
	# e.g., GWD (Gambian) individuals are indicated in the 1KG population
	# description file but are not found in the phenotype file.
	if (sum(fam_tab$V6 == 2) == 0) {
		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 										"Correct spelling or try another population."), pop_two))
	}
	if (sum(fam_tab$V6 == 1) == 0) {
 		 stop(sprintf(paste("No '%s' individuals found in PLINK phenotype file.",
		 										"Correct spelling or try another population."), pop_one))
  }

	# Write out new fam file
	fam_folder <- dirname(genotype_file)
  fam_out <- sprintf("%s/%s.fam", fam_folder, output_file)
  write.table(fam_tab, file=fam_out, row=FALSE, quote=FALSE, col=FALSE)

	n_pop_one <- length(which(fam_tab$V6 == 1))
	n_pop_two <- length(which(fam_tab$V6 == 2))

	cat(sprintf("** Found %s %s and %s %s genotypes\n", n_pop_one, pop_one, n_pop_two, pop_two))
	cat(sprintf("* New population-coded fam file written to %s\n", fam_out))
}
