#' Generates SNP lists per selection-enriched and unenriched pathway
#' as determined by GSEA
#'
#' @param genotype_file (char) path to file with SNP genotype data (PLINK format).
#' @param results_file (char) path to files with GSEA results.
#' 		Strutured to compare results of two population analyses
#' 		i.e., CEU vs. YRI and CEU vs. LWK.
#' @param gseaStat_file (char) path to GSEA statistics file.
#' @param snp2gene_file (char) path to SNP-gene mapping file.
#' @param famF (char) path to PLINK population coded fam file.
#' @param ENRICH_NES (integer) NES cutoff to select validated selection-enriched
#'		pathways (default=3)
#' @param UNENRICH_NES (integer) NES cutoff to select unenriched pathways
#'    (default=0.1)
#' @param enrich_folder (char) path to directory to store output files
#' 		(PLINK files per selection-enriched gene set)
#' @param enrichEM_folder (char) path to directory to store output files
#' 		(PLINK files per selection-enriched pathway, grouped via AutoAnnotate)
#' @param unenrich_folder (char) path to directory to store output files
#' 		(PLINK files per unenriched gene set)
#'
#' @return none
#' @export
#'

writePathFiles <- function(genotype_file, results_file, gseaStat_file,
													 snp2gene_file, famF,
												 	 ENRICH_NES, UNENRICH_NES,
											   	 enrich_folder, enrichEM_folder, unenrich_folder) {
	# Merge GSEA results from both population comparisons and change column names
	# so the results from each comparison can be identified after merging
	pathRes <- lapply(results_file, read.delim)
	for (i in seq_along(pathRes)) {
		colnames(pathRes[[i]])[2:7] <- paste(colnames(pathRes[[i]][, c(2:7)]), i, sep="_")
	}
	res_merge <- join_all(pathRes, by="Geneset")

	pathStats <- function(pathSet, NEScut, output_folder) {
		# Get selection-enriched and unenriched pathways from GSEA results
		cat(sprintf("*Determining %s pathways from GSEA results...\n", pathSet))
		if (pathSet == "enrich" | pathSet == "enrich_EM") {
			cat(sprintf("  Selection-enriched NES threshold: %g\n", NEScut))
			paths <- filter(res_merge, NES_1 >= NEScut & FDR_1 <= 0.05)
			if (length(pathRes) > 1) { # filter by second population analysis if run
				paths <- filter(paths, NES_2 >= NEScut & FDR_2 <= 0.05)
			}
		} else if (pathSet == "unenrich") {
			cat(sprintf("\tUnenriched NES threshold: %g\n", NEScut))
			paths <- filter(res_merge, NES_1 <= NEScut & NES_1 >= -NEScut)
			if (length(pathRes) > 1) { # filter by second population analysis if run
				paths <- filter(paths, NES_2 <= NEScut & NES_2 >= -NEScut)
			}
		}
		write.table(paths, file=sprintf("%s/results_%s.txt", output_folder, pathSet),
								col=TRUE, row=FALSE, quote=FALSE, sep="\t")

		cat(sprintf("\n*Pulling stats for %i pathways...\n", nrow(paths)))
		# Print list of pathway names, subset GSEA stat file with those
		# pathway names and read back into R
		path_names <- paths[,1]
		write.table(path_names,
								file=sprintf("%s/pathways_%s.txt", output_folder, pathSet),
								col=FALSE, row=FALSE, quote=FALSE)

		pathsF <- sprintf("%s/pathways_%s.txt", output_folder, pathSet)

		# Subset gseaStatFile.txt by enriched and unenriched pathways
		## --colour-never flag specifically for OSX, issue with readLines reading
		## ANSI color codes from grep output
		cat(sprintf("*Subsetting GSEA statistics file for %s pathways\n", pathSet))
		system(sprintf("grep --colour=never -f %s %s > %s/gseaStatFile_%s.txt",
							      pathsF, gseaStat_file, output_folder, pathSet))
		statsF <- sprintf("%s/gseaStatFile_%s.txt", output_folder, pathSet)

		no_col <- max(count.fields(statsF)) #get max number of columns in file
		stats <- readLines(statsF)
		stats <- str_split_fixed(stats, "\t", no_col)
		snp_stats <- t(stats) #transpose data

		if (pathSet == "enrich_EM") {
			# Load pathway groupings (written by plotEmap.R)
			load(sprintf("%s/pathway_groups.rda", output_folder))
			for (i in seq_along(pathway_groups)) {
				path <- snp_stats[,which(snp_stats[1,] %in% pathway_groups[[i]])]
				# Remove pathway names
				if (!length(ncol(blah))) {
					path <- path[-1]
				} else {
					path <- path[-1,]
				}
				# Split original df into 3 columns (gene, snp, fst value)
				path <- as.data.frame(str_split_fixed(path, ",", 3))
				# Recode blank cells to NA then remove them
				path[path == ""] <- NA
				path <- na.omit(path)
				# Get unique SNPs/genes in pathway group
				path <- unique(path)
				# Replace spaces with underscore in pathway name strings
				name <- names(pathway_groups)[i]
				name <- gsub('([[:punct:]])|\\s+', '_', name)

				# Separate dfs for snps and genes per pathway
				snp_list <- as.data.frame(path[,2])
				gene_list <- as.data.frame(path[,1])

				snp_file <- file.path(sprintf("%s/%s.snps", output_folder, name))
				cat(sprintf("\n*Generating list for SNPs in %s pathway...", name))
				write.table(snp_list, file=snp_file, col=FALSE, row=FALSE, quote=FALSE)
				cat(" done.\n")

				gene_file <- file.path(sprintf("%s/%s.genes", output_folder, name))
				cat(sprintf("*Generating list for genes in %s pathway...", name))
				write.table(gene_list, file=gene_file, col=FALSE, row=FALSE, quote=FALSE)
				cat(" done.\n")

				outF <- substr(snp_file, 0, nchar(snp_file)-5)
				str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s --extract %s",
												genotype_file, genotype_file, famF, snp_file)
				str2 <- sprintf("--make-bed --allow-no-sex --out %s", outF)
				cmd <- sprintf("%s %s", str1, str2)
				system(cmd)
			}
		} else {
			for (i in 1:ncol(snp_stats)) {
				# Remove pathway names
				path <- snp_stats[-1,i]
				# Split original df into 3 columns (gene, snp, fst value)
				path <- as.data.frame(str_split_fixed(path, ",", 3))
				# Recode blank cells to NA then remove them
				path[path == ""] <- NA
				path <- na.omit(path)

				# Separate dfs for snps and genes per pathway
				snp_list <- as.data.frame(path[,2])
				gene_list <- as.data.frame(path[,1])

				# Replace spaces with underscore in pathway name strings
				name <- gsub(" ", "_", snp_stats[1,i], fixed=TRUE)
				# Replace all special chars with underscore for PLINK compatibility
				name <- gsub('([[:punct:]])|\\s+', '_', name)

				snp_file <- file.path(sprintf("%s/%s.snps", output_folder, name))
				cat(sprintf("\n*Generating list for SNPs in %s pathway...", name))
				write.table(snp_list, file=snp_file, col=FALSE, row=FALSE, quote=FALSE)
				cat(" done.\n")

				gene_file <- file.path(sprintf("%s/%s.genes", output_folder, name))
				cat(sprintf("*Generating list for genes in %s pathway...", name))
				write.table(gene_list, file=gene_file, col=FALSE, row=FALSE, quote=FALSE)
				cat(" done.\n")

				# Subset original PLINK file with each pathway SNP file and write out
				# new set of files per pathway (for use in SNP association calculation)
				outF <- substr(snp_file, 0, nchar(snp_file)-5)
				str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s --extract %s",
												genotype_file, genotype_file, famF, snp_file)
				str2 <- sprintf("--make-bed --allow-no-sex --out %s", outF)
				cmd <- sprintf("%s %s", str1, str2)
				system(cmd)
			}
		}

		# Concatenate SNP files into master file
		concatFiles <- function(x) {
			cat(sprintf("*Combining all %s into master file...", x))
			outF <- sprintf("%s/%s_%s", output_folder, x, pathSet)
			system(sprintf("cat %s/*.%s > %s.txt", output_folder, x, outF))
			cat(sprintf(" file written to %s.txt.\n", outF))

			cat(sprintf("*Writing file with only unique %s...", x))
			dat <- read.table(sprintf("%s.txt", outF), h=FALSE, as.is=TRUE)
			dat_unique <- as.data.frame(unique(dat))
			write.table(dat_unique, file=sprintf("%s_unique.txt", outF),
									col=FALSE, row=FALSE, quote=FALSE)
			cat(sprintf(" file written to %s_unique.txt.\n", outF))
		}
		mapply(concatFiles, c("snps", "genes"))
	}

 # Run function for enriched and unenriched pathways
 set_list <- c("enrich", "enrich_EM", "unenrich")
 nes_list <- c(ENRICH_NES, ENRICH_NES, UNENRICH_NES)
 dir_list <- c(enrich_folder, enrichEM_folder, unenrich_folder)

 mapply(pathStats, set_list, nes_list, dir_list)
}
