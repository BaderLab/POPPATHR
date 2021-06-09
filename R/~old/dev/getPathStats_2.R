#' Generates SNP lists per high-confidence pathway as determined by GSEA

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param resF (char) path to files with GSEA results
#' 		Strutured to compare same population analysis using two datasets
#' 		i.e., CEU vs. ASW using 1KG dataset and HM3 dataset
#' @param gseaStatF (char) path to GSEA statistics file
#' 		(generated via setupGSEArun.R)
#' @param noFilter (char) get all possible genes in a pathway; no MAF filter
#' @param snp2gene (char) path to file with SNP-gene mappings
#' @param NEScut (integer) indicates NES cutoff to filter pathawys
#' @param hcFDRcut_res1 (integer) indicates FDR cutoff to filter significant
#'		high-confidence pathways from the first dataset (e.g., HM3)
#' @param hcFDRcut_res2 (integer) indicates FDR cutoff to filter significant
#'		high-confidence pathways from the second dataset (e.g., PNC)
#' @param PLINK (char) path to PLINK executable
#' @param hcOutDir (char) path to directory to store output files
#' 		(PLINK binary files for each high-confidence pathway)
#' @param lcOutDir (char) path to directory to store output files
#' 		(PLINK binary files for each low-confidence pathway)
#' @return
#' @export

getPathStats <- function(genoF, resF, gseaStatF, PLINK,
									 			 noFilter, snp2gene,
												 NEScut=1.4, FDRcut_res1=0.1, FDRcut_res2=0.05,
											   hcOutDir, lcOutDir) {

	# Merge GSEA results from both datasets and change column names
	# so the results from each dataset can be identified after merging
	pathRes <- lapply(resF, read.delim)
	for (i in 1:length(pathRes)) {
		colnames(pathRes[[i]])[2:7] <- paste(colnames(pathRes[[i]][, c(2:7)]),
																				 sprintf("res%i", i), sep="_") }
	  res_merge <- merge(pathRes[[1]], pathRes[[2]], by="Geneset")

	# Get the high confidence and low confidence pathways from GSEA results
	pathStats <- function(pathSet, outDir) {

		if (noFilter == TRUE) {
			outDir <- sprintf("%s_noMAF", outDir)
			if (!file.exists(outDir)) dir.create(outDir)

			# Get all possible SNP-gene pairings per pathway
			snp2gene <- read.table(snp2geneF, h=F, as.is=T,
														 col.names=c("snp", "gene", "bp")) }

		if (pathSet == "hc") {
			paths <- filter(res_merge,
											NES_res1 >= NEScut & NES_res2 >= NEScut &
											FDR_res1 <= FDRcut_res1 & FDR_res2 <= FDRcut_res2)
			write.table(paths, file=sprintf("%s/pathways_%s.txt", outDir, pathSet),
								  col=T, row=F, quote=F, sep="\t")

		} else if (pathSet == "lc") {
			tmp <- filter(res_merge,
										NES_res1 <= NEScut & NES_res1 >= -NEScut &
										NES_res2 <= NEScut & NES_res2 >= -NEScut)
			tmp <- tmp[order(abs(tmp$NES_res1), abs(tmp$NES_res2)), ]
			hc <- read.delim(sprintf("%s/pathways_hc.txt", hcOutDir), h=T, as.is=T)
			paths <- data.frame()
			paths <- tmp[1,]
			for (i in 1:nrow(tmp)) {
			  if (sum(paths$Size_res1) <= sum(hc$Size_res1)) {
			    paths <- rbind(paths, tmp[i+1,])
			  } else {
			  return()
			  }
			}

			write.table(paths, file=sprintf("%s/pathways_%s.txt", outDir, pathSet),
								  col=T, row=F, quote=F, sep="\t")

		}

		cat(sprintf("\n*Pulling stats for %i pathways...\n", nrow(paths)))

		# Print list of pathway names, subset GSEA stat file with those
		# pathway names and read back into R
		path_names <- paths[,1]
		write.table(path_names,
								file=sprintf("%s/pathway_names_%s.txt", outDir, pathSet),
								col=F, row=F, quote=F)

		pathsF <- sprintf("%s/pathway_names_%s.txt", outDir, pathSet)
		system(sprintf("grep -f %s %s > %s/gseaStatFile_%s.txt",
							      pathsF, gseaStatF, outDir, pathSet))
		statsF <- sprintf("%s/gseaStatFile_%s.txt", outDir, pathSet)

		no_col <- max(count.fields(statsF)) #get max number of columns in stat file
		stats <- readLines(statsF)
		stats <- str_split_fixed(stats, "\t", no_col)
		snp_stats <- t(stats) #transpose data

		for (i in 1:ncol(snp_stats)) {
			#remove first row from each df (pathway name)
			path <- snp_stats[-1,i]
			#split into 3 columns by comma delimiter
			path <- as.data.frame(str_split_fixed(path, ",", 3))
			#recode blank cells to NA then remove them
			path[path == ""] <- NA
			path <- na.omit(path)

			#print out dataframe of each leading edge gene per pathway
			###LE_genes <- as.data.frame(path[1,1])
			###write.table(LE_genes, "LE_genes.txt", col=F, row=F, quote=F, append=T)
			#to remove LE genes from gmt file
			#system(grep -Fvf LE_genes [gmt] > [gmt]_sans_LE)

			if (noFilter == TRUE) {
				#keep column containing genes from each pathway list and merge with
				#snp2gene file (all possible SNP-gene pairings)
				gene_list <- as.data.frame(path[,1])
				colnames(gene_list) <- "gene"
				#now we have all possible SNPs mapped to a single gene
				path <- merge(gene_list, snp2gene, by="gene")
			}

			#replace spaces with underscore in pathway name strings
			names <- gsub(" ", "_", snp_stats[1,i], fixed=T)
			#remove paranthesis and backwards slash from pathway name strings
			#(PLINK doesn't handle them)
			names <- gsub("\\(", "", names)
			names <- gsub("\\)", "", names)
			names <- gsub("\\/", "", names)
			names <- gsub("\\'", "", names) #also remove apostrophes
			names <- gsub("\\,", "", names) #and commas

			#separate dfs for snps and genes per pathway list
			snp_list <- as.data.frame(path[,2])
			gene_list <- as.data.frame(path[,1])

			snp_file <- file.path(sprintf("%s/%s.snps", outDir, names))
			cat(sprintf("\n*Generating list for all SNPs in %s pathway...", names))
			write.table(snp_list, file=snp_file, col=F, row=F, quote=F)
			cat(" done.\n")

			gene_file <- file.path(sprintf("%s/%s.genes", outDir, names))
			cat(sprintf("*Generating list for all genes in %s pathway...", names))
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
	 # Run function for high and low confidence pathways
	 pathStats(pathSet="hc", outDir=hcOutDir)
	 pathStats(pathSet="lc", outDir=lcOutDir)

}
