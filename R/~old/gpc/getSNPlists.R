#' Generates SNP lists per high-confidence pathway as determined by GSEA

#' @param genoF (char) path to file with SNP genotype data (PLINK format)
#' @param resFiles (char) path to files with GSEA results
#' Strutured to compare same population analysis using two datasets
#' i.e., CEU vs. ASW using 1KG dataset and HM3 dataset
#' @param FDRcut (integer) indicates FDR cutoff to filter significant pathways
#' @param PLINK (char) path to PLINK executable
#' @param gseaStatF (char) path to GSEA statistics file (generated via setupGSEArun.R)
#' @param testStatF (char) path to SNP test statistic file (generated via calcMAFdiff.R)
#' @param outDir (char) path to directory to store output files (PLINK binary files for
#' each high-confidence pathway)
#' @export

getSNPlists <- function(genoF, resFiles, FDRcut=0.15,
											  PLINK, gseaStatF, testStatF, outDir) {

	# Merge GSEA results for both datasets
	pathRes <- lapply(resFiles, read.delim)
	for (i in 1:length(pathRes)) {
		colnames(pathRes[[i]])[2:7] <- paste(colnames(pathRes[[i]][, c(2:7)]),
																				 sprintf("res%i", i), sep="_")

	}

	res_merge <- merge(pathRes[[1]], pathRes[[2]], by="Geneset")

	# Get the high confidence pathways from GSEA results
	high_conf <- filter(res_merge, FDR_res1 < FDRcut & FDR_res2 < FDRcut)
	# Print table of high confidence pathway GSEA stats
	write.xlsx(high_conf,
						 file=sprintf("%s/high_conf_pathways.xlsx", outDir),
						 col=T, row=F)

  # Temporarily print list of pathway names for matching
	high_path_names <- high_conf[,1]
  write.table(high_path_names,
							file=sprintf("%s/hc_genesets.tmp", outDir),
							col=F, row=F, quote=F)

	# Subset GSEA stat file with high-conf pathway names and read back into R
	hcPathsF <- sprintf("%s/hc_genesets.tmp", outDir)
	system(sprintf("grep -f %s %s > %s/hc_stats.tmp",
						      hcPathsF, gseaStatF, outDir))
	hcStatsF <- sprintf("%s/hc_stats.tmp", outDir)
	no_col <- max(count.fields(hcStatsF)) #get max number of columns in stat file

	# PROBLEM WITH THIS STEP: for some reason when reading the stat file back
	# into R, it's cutting some pathways out... tried read.delim() and different
	# delimiter but doesn't work...
	hcStats <- read.table(hcStatsF, sep="\t", fill=T, col.names=1:no_col)
	hc_snp_stats <- t(hcStats) #transpose data

	for (i in 1:ncol(hc_snp_stats)) {
		#remove hc row from each df; the pathway name - re-formatted later
		hc_path <- hc_snp_stats[-1,i]
		#split into 3 columns by comma delimiter
		hc_path <- as.data.frame(str_split_fixed(hc_path, ",", 3))
		#recode blank cells to NA then remove them
		hc_path[hc_path == ""] <- NA
		hc_path <- na.omit(hc_path)
		#keep column containing rsIDs from each pathway list
		hc_snp_list <- as.data.frame(hc_path[,2])
		#replace spaces with underscore in pathway name strings
		hc_names <- gsub(" ", "_", hc_snp_stats[1,i], fixed=T)
		#remove paranthesis and backwards slash from pathway name strings
		#(PLINK doesn't handle them)
		hc_names <- gsub("\\(", "", hc_names)
		hc_names <- gsub("\\)", "", hc_names)
		hc_names <- gsub("\\/", "", hc_names)

		path_file <- file.path(sprintf("%s/%s.snps", outDir, hc_names))
		cat(sprintf("Generating SNP list for %s pathway...", hc_names))
		write.table(hc_snp_list, file=path_file, col=F, row=F, quote=F)
		cat(" done.\n")

		# Subset original PLINK file with each pathway SNP file
		str1 <- sprintf("%s --bfile %s --extract %s", PLINK, genoF, path_file)
		str2 <- sprintf("--make-bed --allow-no-sex --out %s",
										file_path_sans_ext(path_file))
		cmd <- sprintf("%s %s", str1, str2)
		system(cmd)
	}

	# Concatenate SNP files into master file
	cat("Combining all SNPs into master SNP file...")
	system(sprintf("cat %s/*.snps > %s/all_snps.txt", outDir, outDir))
	cat(sprintf(" file written to %s/all_snps.txt.\n", outDir))

	snps <- read.table(sprintf("%s/all_snps.txt", outDir), h=F, as.is=T)
	snps_unique <- as.data.frame(unique(snps))
	cat("Writing file with only unique SNPs...")
	write.table(snps_unique,
							file=sprintf("%s/all_snps_unique.txt", outDir),
							col=F, row=F, quote=F)
	cat(sprintf(" file written to %s/all_snps_unique.txt.\n", outDir))

	cat("\n---------------------------------------------------------------------")
	cat(sprintf(paste("\n%i total SNPs (%i unique SNPs) mapped to genes",
										"within\n%i high confidence pathways at FDR < %g.\n"),
									 	nrow(snps), nrow(unique(snps)), ncol(hc_snp_stats), FDRcut))
	cat("---------------------------------------------------------------------\n")

	# Clean up
	unlink(c(hcPathsF, hcStatsF))
}
