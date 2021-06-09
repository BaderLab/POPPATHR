#' Generates SNP lists per high-confidence pathway as determined by GSEA

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param resF (char) path to files with GSEA results
#' 		Strutured to compare same population analysis using two datasets
#' 		i.e., CEU vs. ASW using 1KG dataset and HM3 dataset
#' @param gseaStatF (char) path to GSEA statistics file
#' 		(generated via setupGSEArun.R)
#' @param snp2gene (char) path to SNP-gene mapping file
#' @param LEsubset (logical) only include leading edge genes in pathway files
#' 		(default=FALSE)
#' @param hcFDRcut (integer) indicates FDR cutoff to select significant
#'		high-confidence pathways from both datasets (default=0.1)
#' @param lcNEScut (integer) indicates NES cutoff to select low-confidence
#'    pathways from both datasets (default=0.1)
#' @param PLINK (char) path to PLINK executable
#' @param hcOutDir (char) path to directory to store output files
#' 		(PLINK binary files for each high-confidence pathway)
#' @param lcOutDir (char) path to directory to store output files
#' 		(PLINK binary files for each low-confidence pathway)
#' @return
#' @export

getPathStats <- function(genoF, resF, gseaStatF, snp2gene, PLINK,
												 LEsubset=FALSE, hcFDRcut=0.1, lcNEScut=0.1,
											   hcOutDir, lcOutDir) {

	# Merge GSEA results from both datasets and change column names
	# so the results from each dataset can be identified after merging
	pathRes <- lapply(resF, read.delim)
	for (i in 1:length(pathRes)) {
		colnames(pathRes[[i]])[2:7] <- paste(colnames(pathRes[[i]][, c(2:7)]),
																				 sprintf("res%i", i), sep="_") }

	res_merge <- merge(pathRes[[1]], pathRes[[2]], by="Geneset")

	# Get the high confidence and low confidence pathways from GSEA results
	pathStats <- function(pathSet, FDRcut, NEScut, outDir) {

		if (pathSet == "hc") {
			cat(sprintf("*Determining %s pathways from GSEA results...\n", pathSet))
			cat(sprintf("\tHigh-confidence FDR cutoff: %g\n", FDRcut))

			paths <- filter(res_merge,
											FDR_res1 <= FDRcut & FDR_res2 <= FDRcut)
			write.table(paths, file=sprintf("%s/results_%s.txt", outDir, pathSet),
								  col=T, row=F, quote=F, sep="\t")

		} else if (pathSet == "lc") {
			cat(sprintf("*Determining %s pathways from GSEA results...\n", pathSet))
			cat(sprintf("\tLow-confidence NES cutoff: %g\n", NEScut))

			paths <- filter(res_merge,						#filter by NES
											NES_res1 <= NEScut & NES_res1 >= -NEScut &
											NES_res2 <= NEScut & NES_res2 >= -NEScut)
			write.table(paths, file=sprintf("%s/results_%s.txt", outDir, pathSet),
								  col=T, row=F, quote=F, sep="\t")

		} else if (pathSet == "all") {
			paths <- pathRes[[1]]
		}

		cat(sprintf("\n*Pulling stats for %i %s pathways...\n",
				nrow(paths), pathSet))

		# Print list of pathway names, subset GSEA stat file with those
		# pathway names and read back into R
		path_names <- paths[,1]
		write.table(path_names,
								file=sprintf("%s/pathways_%s.txt", outDir, pathSet),
								col=F, row=F, quote=F)

		pathsF <- sprintf("%s/pathways_%s.txt", outDir, pathSet)

		# Subset gseaStatFile.txt by hc and lc pathways
		## --colour-never flag specifically for OSX, issue with readLines reading
		## ANSI color codes from grep output
		system(sprintf("grep --colour=never -f %s %s > %s/gseaStatFile_%s.txt",
							      pathsF, gseaStatF, outDir, pathSet))
		statsF <- sprintf("%s/gseaStatFile_%s.txt", outDir, pathSet)

		if (LEsubset == TRUE) {
			cat("*OPTIONAL: subsetting pathways by leading edge genes...\n")
			# Create new output directory to list pathway files with only LE genes
			outDir <- sprintf("%s/%s_snps_LE", ldDir, pathSet)
			if (!file.exists(outDir)) dir.create(outDir)

			# Subset gseaLEout.txt by hc and lc pathways
			system(sprintf("grep -f %s %s > %s/gseaLEout_%s.txt",
											pathsF, gseaLEoutF, outDir, pathSet))
			LEoutF <- sprintf("%s/gseaLEout_%s.txt", outDir, pathSet)

			# read leading edge genes
		  no_col <- max(count.fields(LEoutF)) #get max number of columns in file
		  le <- readLines(LEoutF)
		  le <- str_split_fixed(le, "\t", no_col)
		  le_stats <- t(le) #transpose data
		  le_stats <- le_stats[-c(1,2),] #remove pathway names
		}

		no_col <- max(count.fields(statsF)) #get max number of columns in file
		stats <- readLines(statsF)
		stats <- str_split_fixed(stats, "\t", no_col)
		snp_stats <- t(stats) #transpose data

		for (i in 1:ncol(snp_stats)) {

			if (LEsubset == TRUE) {
				# Remove pathway names
				path <- snp_stats[-1,]

				le_genes <- le_stats[,i]
				le_genes[le_genes == ""] <- NA
				le_genes <- na.omit(le_genes)

				# Paste all LE genes together sep by "|" in order to subset those
				# genes from the original gseaStatFile.txt
				le_genes <- paste(le_genes, collapse="|")
				le_sub <- grep(le_genes, path[,i], value=T)

				# Split df into 3 columns by comma delimiter
				le_split <- as.data.frame(str_split_fixed(le_sub, ",", 3))

				# Separate dfs for snps and genes per pathway list
				snp_list <- as.data.frame(le_split[,2])
				gene_list <- as.data.frame(le_split[,1])
			} else {
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
			}

			# Replace spaces with underscore in pathway name strings
			names <- gsub(" ", "_", snp_stats[1,i], fixed=T)
			# Replace all special chars with underscore for PLINK compatibility
			names <- gsub('([[:punct:]])|\\s+', '_', names)

			snp_file <- file.path(sprintf("%s/%s.snps", outDir, names))
			cat(sprintf("\n*Generating list for SNPs in %s pathway...", names))
			write.table(snp_list, file=snp_file, col=F, row=F, quote=F)
			cat(" done.\n")

			gene_file <- file.path(sprintf("%s/%s.genes", outDir, names))
			cat(sprintf("*Generating list for genes in %s pathway...", names))
			write.table(gene_list, file=gene_file, col=F, row=F, quote=F)
			cat(" done.\n")

			# Subset original PLINK file with each pathway SNP file and write out
			# new sets of PLINK files per pathway
			str1 <- sprintf("%s --bed %s.bed --bim %s.bim --fam %s --extract %s",
											PLINK, genoF, genoF, realFAM, snp_file)
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
		snps <- read.table(sprintf("%s/snps_%s.txt", outDir, pathSet), h=F, as.is=T)
		snps_unique <- as.data.frame(unique(snps))
		write.table(snps_unique,
								file=sprintf("%s/snps_unique_%s.txt", outDir, pathSet),
								col=F, row=F, quote=F)
		cat(sprintf(" file written to %s/snps_unique_%s.txt.\n", outDir, pathSet))
	}

 # Run function for high-confidence pathways
 pathStats(pathSet="hc",
	 				 FDRcut=hcFDRcut,
					 outDir=hcOutDir)

 # Run function for low-confidence pathways
 pathStats(pathSet="lc",
 					 NEScut=lcNEScut,
 					 outDir=lcOutDir)

}
