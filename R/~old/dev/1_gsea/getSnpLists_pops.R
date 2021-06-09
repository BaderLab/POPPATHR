#-----------------------------INPUT FILES---------------------------------------
dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways"

plinkPop1F <- sprintf(paste0("%s/data/PNC/CEU/PNC_imputed_merged.CLEAN_FINAL",
					     						   "_5sd_noAxiom_CEU"), dataDir)
plinkPop2F <- sprintf(paste0("%s/data/PNC/ASW/PNC_imputed_merged.CLEAN_FINAL",
					     						   "_5sd_noAxiom_ASW"), dataDir)

gseaDir1 <- sprintf(paste0("%s/methods/1_gsea/PNC/CEU_ASW/",
													 "out_161109_MAFdiff_CEU-ASW"), dataDir)
gseaDir2 <- sprintf(paste0("%s/methods/1_gsea/HM3/CEU_ASW/",
													 "out_161109_MAFdiff_CEU-ASW"), dataDir)

pathwayRes1 <- sprintf("%s/gsea/results.txt", gseaDir1)
pathwayRes2 <- sprintf("%s/gsea/results.txt", gseaDir2)
gseaStatF <- sprintf("%s/gsea/gseaStatFile.txt", gseaDir1)
testStatF <- sprintf("%s/gsea/markerMAF.txt", gseaDir1)

outDir  <- sprintf("%s/input_snps", dataDir)
		highPop1Dir <- sprintf("%s/high_conf/pnc/ceu", outDir)
		highPop2Dir <- sprintf("%s/high_conf/pnc/asw", outDir)
		randPop1Dir <- sprintf("%s/rand/pnc/ceu", outDir)
		randPop2Dir <- sprintf("%s/rand/pnc/asw", outDir)

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
if (!file.exists(highPop1Dir)) dir.create(highPop1Dir)
if (!file.exists(highPop2Dir)) dir.create(highPop2Dir)
if (!file.exists(randpop1Dir)) dir.create(randpop1Dir)
if (!file.exists(randpop2Dir)) dir.create(randpop2Dir)

#------------------------------------------------
getSnpList <- function(plinkPop1F, plinkPop2F,
											 pathwayRes1, pathwayRes2,
											 gseaStatF, testStatF) {

	# Merge GSEA results for both datasets
  res1 <- read.delim(pathwayRes1, h=T, as.is=T)
  res2 <- read.delim(pathwayRes2, h=T, as.is=T)
  colnames(res1)[2:7] <- paste(colnames(res1[,c(2:7)]), "res1", sep="_")
  colnames(res2)[2:7] <- paste(colnames(res2[,c(2:7)]), "res2", sep="_")
  res.merge <- merge(res1, res2, by="Geneset")

	# Get the high confidence pathways from GSEA results
	hc.res <- which(res.merge$FDR_res1 < 0.05 & res.merge$FDR_res2 < 0.1)
  hc.res <- res.merge[hc.res, ]
	# Print table of high confidence pathway GSEA stats
	write.xlsx(hc.res,
						 file=sprintf("%s/high_conf_pathways.xlsx", outDir),
						 col=T, row=F)
  hc.path.names <- hc.res[,1]
	# Temporarily print list of pathway names for matching
  write.table(hc.path.names,
							file=sprintf("%s/hc_genesets.tmp", outDir),
							col=F, row=F, quote=F)

	# Subset GSEA stat file with high-conf pathway names and read back into R
	hcPathsF <- sprintf("%s/hc_genesets.tmp", outDir)
	system(sprintf("grep -f %s %s > %s/hc_stats.tmp",
						      hcPathsF, gseaStatF, outDir))
	hcStatsF <- sprintf("%s/hc_stats.tmp", outDir)
	no.col <- max(count.fields(hcStatsF, sep="\t"))
  hc.stats <- read.table(hcStatsF, sep="\t", fill=T, col.names=1:no.col)
	hc.snp.stats <- t(hc.stats) #transpose data

	for (i in 1:ncol(hc.snp.stats)) {
		#remove top row from each df; the pathway name - re-formatted later
		hc.path <- hc.snp.stats[-1,i]
		#split into 3 columns by comma delimiter
		hc.path <- as.data.frame(str_split_fixed(hc.path, ",", 3))
		#recode blank cells to NA then remove them
		hc.path[hc.path == ""] <- NA
		hc.path <- na.omit(hc.path)
		#keep column containing rsIDs from each pathway list
		hc.snp.list <- as.data.frame(hc.path[,2])
		#replace spaces with underscore in pathway name strings
		hc.names <- gsub(" ", "_", hc.snp.stats[1,i], fixed=T)
		#remove paranthesis from pathway name strings (PLINK doesn't handle them)
		hc.names <- gsub("\\(", "", hc.names)
		hc.names <- gsub("\\)", "", hc.names)

		cat(sprintf("Generating SNP list for %s pathway...", hc.names))

		snps.pop1 <- file.path(sprintf("%s/%s.snps", highPop1Dir, hc.names))
		write.table(hc.snp.list, file=snps.pop1, col=F, row=F, quote=F)

		snps.pop2 <- file.path(sprintf("%s/%s.snps", highPop2Dir, hc.names))
		write.table(hc.snp.list, file=snps.pop2, col=F, row=F, quote=F)

		cat(" done.\n")

		# Subset original PLINK file with each pathway SNP file
		## Population 1 -- CEU
		cat(sprintf("PLINK subset population 1 (%s)\n", basename(plinkPop1F)))
		str1 <- sprintf("%s --bfile %s --extract %s", PLINK, plinkPop1F, snps.pop1)
		str2 <- sprintf("--make-bed --allow-no-sex --out %s",
										file_path_sans_ext(snps.pop1))
		cmd <- sprintf("%s %s", str1, str2)
		system(cmd)

		# Subset original PLINK file with each pathway SNP file
		## Population 2 -- ASW
		cat(sprintf("PLINK subset population 1 (%s)\n", basename(plinkPop1F)))
		str1 <- sprintf("%s --bfile %s --extract %s", PLINK, plinkPop2F, snps.pop2)
		str2 <- sprintf("--make-bed --allow-no-sex --out %s",
										file_path_sans_ext(snps.pop2))
		cmd <- sprintf("%s %s", str1, str2)
		system(cmd)
	}

	# Concatenate SNP files into master file
	system(sprintf("cat %s/*.snps > %s/all_snps.txt", highPop1Dir, highPop1Dir))
	system(sprintf("cat %s/*.snps > %s/all_snps.txt", highPop2Dir, highPop2Dir))

	all.snps.pop1 <- read.table(sprintf("%s/all_snps.txt", highPop1Dir),
																			h=F, as.is=T)
	all.snps.pop2 <- read.table(sprintf("%s/all_snps.txt", highPop2Dir),
																			h=F, as.is=T)

	# Get only unique SNPs
	snps.unique.pop1 <- as.data.frame(unique(all.snps.pop1))
	snps.unique.pop2 <- as.data.frame(unique(all.snps.pop2))

	write.table(snps.unique.pop1,
							file=sprintf("%s/all_snps_unique.txt", highPop1Dir),
							col=F, row=F, quote=F)

	write.table(snps.unique.pop2,
							file=sprintf("%s/all_snps_unique.txt", highPop2Dir),
							col=F, row=F, quote=F)

	# Randomly take sample SNPs from original PLINK genotype file
	## for use in downstream analyses (using UNIX 'shuf' command)
	rand.pop1 <- sprintf("%s/random_sample_ceu.snps", randpop1Dir)
	rand.pop2 <- sprintf("%s/random_sample_asw.snps", randpop2Dir)

	# Take same amount of unique SNPs relative to high conf pathways
	cat(sprintf("\nRandomly sampling %i SNPs from %s...\n",
						  nrow(unique(all.snps.pop1)), basename(plinkPop1F)))
	cat(sprintf("\nRandomly sampling %i SNPs from %s...\n",
						  nrow(unique(all.snps.pop2)), basename(plinkPop2F)))

	# Shuffle PLINK bim file and print 2nd row (rsID)
	system(sprintf("shuf -n %i %s.bim | awk '{print$2}' > %s",
								 nrow(unique(all.snps.pop1)), plinkPop1F, rand.pop1))
	system(sprintf("shuf -n %i %s.bim | awk '{print$2}' > %s",
	  						 nrow(unique(all.snps.pop2)), plinkPop2F, rand.pop2))

	# Generate PLINK binary files for random SNP list
	## Population 1 -- CEU
	str1 <- sprintf("%s --bfile %s --extract %s --make-bed",
									PLINK, plinkPop1F, rand.pop1)
	str2 <- sprintf("--allow-no-sex --out %s", file_path_sans_ext(rand.pop1))
	cmd <- sprintf("%s %s", str1, str2)
	system(cmd)
	cat(" done.\n")

	## Population 2 -- ASW
	str1 <- sprintf("%s --bfile %s --extract %s --make-bed",
									PLINK, plinkPop2F, rand.pop2)
	str2 <- sprintf("--allow-no-sex --out %s", file_path_sans_ext(rand.pop2))
	cmd <- sprintf("%s %s", str1, str2)
	system(cmd)
	cat(" done.\n")

	# Clean up
	unlink(c(hcPathsF, hcStatsF))
}
#------------------------------------------------
getSnpList(plinkPop1F, plinkPop2F,
					 pathwayRes1, pathwayRes2,
					 gseaStatF, testStatF)
