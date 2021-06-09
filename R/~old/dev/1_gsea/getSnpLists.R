#-----------------------------INPUT FILES---------------------------------------
dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways"

plinkF <- sprintf("%s/data/HM3/all/HM3_pops_hg19", dataDir)
#plinkF <- sprintf(paste0("%s/data/PNC/CEU/PNC_imputed_merged.CLEAN_FINAL",
#											   "_5sd_noAxiom_CEU"), dataDir)

gseaDir1 <- sprintf(paste0("%s/methods/1_gsea/PNC/CEU_ASW/",
													 "out_161109_MAFdiff_CEU-ASW"), dataDir)
gseaDir2 <- sprintf(paste0("%s/methods/1_gsea/HM3/CEU_ASW/",
													 "out_161109_MAFdiff_CEU-ASW"), dataDir)

pathwayRes1 <- sprintf("%s/gsea/results.txt", gseaDir1)
pathwayRes2 <- sprintf("%s/gsea/results.txt", gseaDir2)
gseaStatF <- sprintf("%s/gsea/gseaStatFile.txt", gseaDir1)
testStatF <- sprintf("%s/gsea/markerMAF.txt", gseaDir1)

outDir  <- sprintf("%s/input_snps", dataDir)
		highDir <- sprintf("%s/high_conf/hm3", outDir)
		randDir <- sprintf("%s/rand/hm3", outDir)

PLINK <- "/media/catherine/DATAPART1/Software/plink_linux_1.90/plink"

#------------------------------------------------
require(stringr)
require(tools)
require(xlsx)

if (file.exists(outDir)) {
  cat("Directory exists! Not overwriting\n")
  Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)
if (!file.exists(highDir)) dir.create(highDir)
if (!file.exists(randDir)) dir.create(randDir)

#------------------------------------------------
getSnpList <- function(plinkF, pathwayRes1, pathwayRes2,
											 gseaStatF, testStatF) {

	# Merge GSEA results for both datasets
  res1 <- read.delim(pathwayRes1, h=T, as.is=T)
  res2 <- read.delim(pathwayRes2, h=T, as.is=T)
  colnames(res1)[2:7] <- paste(colnames(res1[,c(2:7)]), "res1", sep="_")
  colnames(res2)[2:7] <- paste(colnames(res2[,c(2:7)]), "res2", sep="_")
  res_merge <- merge(res1, res2, by="Geneset")

	# Get the high confidence pathways from GSEA results
	high_conf <- which(res_merge$FDR_res1 < 0.05 & res_merge$FDR_res2 < 0.1)
  high_conf <- res_merge[high_conf, ]
	# Print table of high confidence pathway GSEA stats
	write.xlsx(high_conf,
						 file=sprintf("%s/high_conf_pathways.xlsx", outDir),
						 col=T, row=F)
  high_path_names <- high_conf[,1]
	# Temporarily print list of pathway names for matching
  write.table(high_path_names,
							file=sprintf("%s/top_genesets.tmp", outDir),
							col=F, row=F, quote=F)

	# Subset GSEA stat file with high-conf pathway names and read back into R
	topPathsF <- sprintf("%s/top_genesets.tmp", outDir)
	system(sprintf("grep -f %s %s > %s/top_stats.tmp",
						      topPathsF, gseaStatF, outDir))
	topStatsF <- sprintf("%s/top_stats.tmp", outDir)
	no_col <- max(count.fields(topStatsF, sep="\t"))
  topStats <- read.table(topStatsF, sep="\t", fill=T, col.names=1:no_col)
	top_snp_stats <- t(topStats) #transpose data

	for (i in 1:ncol(top_snp_stats)) {
		#remove top row from each df; the pathway name - re-formatted later
		top_path <- top_snp_stats[-1,i]
		#split into 3 columns by comma delimiter
		top_path <- as.data.frame(str_split_fixed(top_path, ",", 3))
		#recode blank cells to NA then remove them
		top_path[top_path == ""] <- NA
		top_path <- na.omit(top_path)
		#keep column containing rsIDs from each pathway list
		top_snp_list <- as.data.frame(top_path[,2])
		#replace spaces with underscore in pathway name strings
		top_names <- gsub(" ", "_", top_snp_stats[1,i], fixed=T)
		#remove paranthesis from pathway name strings (PLINK doesn't handle them)
		top_names <- gsub("\\(", "", top_names)
		top_names <- gsub("\\)", "", top_names)
		path_file <- file.path(sprintf("%s/%s.snps", highDir, top_names))
		cat(sprintf("Generating SNP list for %s pathway...", top_names))
		write.table(top_snp_list, file=path_file, col=F, row=F, quote=F)
		cat(" done.\n")

		# Subset original PLINK file with each pathway SNP file
		str1 <- sprintf("%s --bfile %s --extract %s", PLINK, plinkF, path_file)
		str2 <- sprintf("--make-bed --allow-no-sex --out %s",
										file_path_sans_ext(path_file))
		cmd <- sprintf("%s %s", str1, str2)
		system(cmd)
	}

	# Concatenate SNP files into master file
	cat("Combining all SNPs into master SNP file...")
	system(sprintf("cat %s/*.snps > %s/all_snps.txt", highDir, highDir))
	cat(sprintf(" file written to %s/all_snps.txt.\n", highDir))

	snps <- read.table(sprintf("%s/all_snps.txt", highDir), h=F, as.is=T)
	snps_unique <- as.data.frame(unique(snps))
	cat("Writing file with only unique SNPs...")
	write.table(snps_unique,
							file=sprintf("%s/all_snps_unique.txt", highDir),
							col=F, row=F, quote=F)
	cat(sprintf(" file written to %s/all_snps_unique.txt.\n", highDir))

	cat("\n---------------------------------------------------------------------")
	cat(sprintf(paste("\n%i total SNPs and %i unique SNPs mapped to genes",
										"within\n%i high confidence pathways.\n"),
									 	nrow(snps), nrow(unique(snps)), ncol(top_snp_stats)))
	cat("---------------------------------------------------------------------\n")
	Sys.sleep(5)

	# Randomly take sample SNPs from original PLINK genotype file
	## for use in downstream analyses (using UNIX 'shuf' command)
	rand_sample <- sprintf("%s/random_sample.snps", randDir)
	# Take same amount of unique SNPs relative to high conf pathways
	cat(sprintf("\nRandomly sampling %i SNPs from %s...\n",
						  nrow(unique(snps)), basename(plinkF)))
	# Shuffle PLINK bim file and print 2nd row (rsID)
  system(sprintf("shuf -n %i %s.bim | awk '{print$2}' > %s",
								 nrow(unique(snps)), plinkF, rand_sample))
	# Generate PLINK binary files for random SNP list
	str1 <- sprintf("%s --bfile %s --extract %s --make-bed",
									PLINK, plinkF, rand_sample)
	str2 <- sprintf("--allow-no-sex --out %s", file_path_sans_ext(rand_sample))
	cmd <- sprintf("%s %s", str1, str2)
	system(cmd)
	cat(" done.\n")

	# Clean up
	unlink(c(topPathsF, topStatsF))
}
#------------------------------------------------
getSnpList(plinkF, pathwayRes1, pathwayRes2,
					 gseaStatF, testStatF)
