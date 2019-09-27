#' Generates SNP lists per selection-enriched and unenriched pathway
#' as determined by GSEA
#'
#' @param genoF (char) path to file with SNP genotype data (PLINK format).
#' @param resF (char) path to files with GSEA results.
#' 		Strutured to compare results of two population analyses
#' 		i.e., CEU vs. YRI and CEU vs. LWK.
#' @param gseaStatF (char) path to GSEA statistics file.
#' @param snp2geneF (char) path to SNP-gene mapping file.
#' @param realFam (char) path to PLINK population coded fam file.
#' @param enrichNES (integer) NES cutoff to select validated selection-enriched
#'		pathways (default=0.3)
#' @param unenrichNES (integer) NES cutoff to select unenriched pathways
#'    (default=0.1)
#' @param enrichDir (char) path to directory to store output files
#' 		(PLINK files per validated selection-enriched pathway)
#' @param unenrichDir (char) path to directory to store output files
#' 		(PLINK files per unenriched pathway)
#'
#' @return none
#' @export
#'

writePathFiles <- function(genoF, resF, gseaStatF, snp2geneF, realFam,
												 	 enrichNES=0.3, unenrichNES=0.1,
											   	 enrichDir, unenrichDir) {
	# Merge GSEA results from both datasets and change column names
	# so the results from each dataset can be identified after merging
	pathRes <- lapply(resF, read.delim)
	for (i in seq_along(pathRes)) {
		colnames(pathRes[[i]])[2:7] <- paste(colnames(pathRes[[i]][, c(2:7)]),
																				 sprintf("res%i", i), sep="_")
	}

	res_merge <- join(pathRes[[1]], pathRes[[2]], by="Geneset")

	pathStats <- function(pathSet, NEScut, outDir) {
		# Get selection-enriched and unenriched pathways from GSEA results
		if (pathSet == "enrich") {
			cat(sprintf("*Determining %sed pathways from GSEA results...\n", pathSet))
			cat(sprintf("  Selection-enriched NES threshold: %g\n", NEScut))
			paths <- filter(res_merge,
											NES_res1 >= NEScut & NES_res2 >= NEScut)
		} else if (pathSet == "unenrich") {
			cat(sprintf("*Determining %sed pathways from GSEA results...\n", pathSet))
			cat(sprintf("\tUnenriched NES threshold: %g\n", NEScut))
			paths <- filter(res_merge,
											NES_res1 <= NEScut & NES_res1 >= -NEScut &
											NES_res2 <= NEScut & NES_res2 >= -NEScut)
		}
		write.table(paths, file=sprintf("%s/results_%s.txt", outDir, pathSet),
								col=TRUE, row=FALSE, quote=FALSE, sep="\t")

		cat(sprintf("\n*Pulling stats for %i %sed pathways...\n", nrow(paths), pathSet))

		# Print list of pathway names, subset GSEA stat file with those
		# pathway names and read back into R
		path_names <- paths[,1]
		write.table(path_names,
								file=sprintf("%s/pathways_%s.txt", outDir, pathSet),
								col=FALSE, row=FALSE, quote=FALSE)

		pathsF <- sprintf("%s/pathways_%s.txt", outDir, pathSet)

		# Subset gseaStatFile.txt by enriched and unenriched pathways
		## --colour-never flag specifically for OSX, issue with readLines reading
		## ANSI color codes from grep output
		system(sprintf("grep --colour=never -f %s %s > %s/gseaStatFile_%s.txt",
							      pathsF, gseaStatF, outDir, pathSet))
		statsF <- sprintf("%s/gseaStatFile_%s.txt", outDir, pathSet)

		no_col <- max(count.fields(statsF)) #get max number of columns in file
		stats <- readLines(statsF)
		stats <- str_split_fixed(stats, "\t", no_col)
		snp_stats <- t(stats) #transpose data

		for (i in 1:ncol(snp_stats)) {
			# Remove pathway names
			path <- snp_stats[-1,i]
			# Split original df into 3 columns by comma delimiter
			path <- as.data.frame(str_split_fixed(path, ",", 3))
			# Recode blank cells to NA then remove them
			path[path == ""] <- NA
			path <- na.omit(path)

			# Separate dfs for snps and genes per pathway list
			snp_list <- as.data.frame(path[,2])
			gene_list <- as.data.frame(path[,1])

			# Replace spaces with underscore in pathway name strings
			names <- gsub(" ", "_", snp_stats[1,i], fixed=TRUE)
			# Replace all special chars with underscore for PLINK compatibility
			names <- gsub('([[:punct:]])|\\s+', '_', names)

			snp_file <- file.path(sprintf("%s/%s.snps", outDir, names))
			cat(sprintf("\n*Generating list for SNPs in %s pathway...", names))
			write.table(snp_list, file=snp_file, col=F, row=F, quote=F)
			cat(" done.\n")

			gene_file <- file.path(sprintf("%s/%s.genes", outDir, names))
			cat(sprintf("*Generating list for genes in %s pathway...", names))
			write.table(gene_list, file=gene_file, col=FALSE, row=FALSE, quote=FALSE)
			cat(" done.\n")

			# Subset original PLINK file with each pathway SNP file and write out
			# new set of files per pathway
			str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s --extract %s",
											genoF, genoF, realFam, snp_file)
			str2 <- sprintf("--make-bed --allow-no-sex --out %s",
											file_path_sans_ext(snp_file))
			cmd <- sprintf("%s %s", str1, str2)
			system(cmd)
		}

		# Concatenate SNP files into master file
		cat("*Combining all SNPs into master SNP file...")
		system(sprintf("cat %s/*.snps > %s/snps_%s.txt", outDir, outDir, pathSet))
		cat(sprintf(" file written to %s/snps_%s.txt.\n", outDir, pathSet))

		cat("*Writing file with only unique SNPs...")
		snps <- read.table(sprintf("%s/snps_%s.txt", outDir, pathSet), h=FALSE, as.is=TRUE)
		snps_unique <- as.data.frame(unique(snps))
		write.table(snps_unique,
								file=sprintf("%s/snps_unique_%s.txt", outDir, pathSet),
								col=FALSE, row=FALSE, quote=FALSE)
		cat(sprintf(" file written to %s/snps_unique_%s.txt.\n", outDir, pathSet))
	}
 # Run function for enriched and unenriched pathways
 set_list <- c("enrich", "unenrich")
 nes_list <- c(enrichNES, unenrichNES)
 dir_list <- c(enrichDir, unenrichDir)

 mapply(pathStats, set_list, nes_list, dir_list)
}
