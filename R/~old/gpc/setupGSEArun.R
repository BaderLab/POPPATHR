#' Sets up and runs GSEA using difference in minor allele frequency (MAF)
#' table containing GSEA results (Geneset, Size, ES, NES, NominalP, FDR, FWER)
#'
#' @param realF (char) path to file with real SNP association statistics
#' @param phenoPerm (logical) runs GSEA using phenotype permutations
#'		(i.e., randomly re-shuffing case/control labels n times).
#' 		if FALSE, will run GSEA using genotype permutations (see setupGSEArun.R)
#' @param permMAF_F (char) path to files with association statistics from
#' 		permuted case/control labels. If phenoPerm is FALSE, leave this blank.
#' @param pathFile (char) path to pathway definitions GMT file
#' @param snp2geneF (char) path to file with snp2gene mappings. Output of
#' 		mapSNP2gene() (found in GWAS2Pathway)
#' @param snp2genedist (integer) value for GSEA --distance.
#'		Max. distance between SNP and gene for their association in snp2geneF
#' @param minGene (integer) value for GSEA --setmin.
#' 		Min. number of genes in a gene set to be considered
#' @param maxGene (integer) value for GSEA --setmax.
#' 		Max. number of genes in a gene set to be considered
#' @param setSeed (integer) value for GSEA --seed
#' @param CALC_GSEA (char) path to calculate_gsea.pl
#' @param COMB_GSEA (char) path to combine_gsea.pl
#' @return (char) list object with the GSEA output filenames
#' @export

require(dplyr)

args       <- commandArgs(TRUE)
realF      <- args[1]
#permMAF_F  <- args[2]
pathFile   <- args[2]
snp2geneF  <- args[3]
CALC_GSEA  <- args[4]
COMB_GSEA  <- args[5]
outDir     <- args[6]

phenoPerm    <- FALSE
snp2genedist <- 10000L
setSeed      <- 42L
minGene      <- 20L
maxGene      <- 200L

# Setup and run GSEA
# NOTE (adapted from the GenGen webpage) Re:importance of perm files
# In the absence of a permutation association results file, the command will
# run a genotype-based permutation, that is, permuting the p-values and
# delta MAF values between SNPs. This is not recommended for SNP arrays
# with high marker-marker linkage disequilibrium patterns, as such
# permutation disrupts the LD structure and may generate biased results.
# Whenever possible, please use raw genotype data to perform phenotype
# permutation and re-calculate SNP association statistics.
phenoPerm <- paste("perl %s %s %s --mapfile %s --permfile %s",
									 "--setstatfile %s --leout %s --setmin %i --setmax %i",
									 "--seed %i --distance %i --log %s/calculate_gsea.log")

genoPerm <- paste("perl %s %s %s --mapfile %s",
									"--setstatfile %s --leout %s --setmin %i --setmax %i",
									"--seed %i --distance %i --log %s/calculate_gsea.log")

statOut <- sprintf("%s/gseaStatFile.txt", outDir)
leOut	  <- sprintf("%s/gseaLEout.txt", outDir)
resOut  <- sprintf("%s/results.txt", outDir)

if (phenoPerm == TRUE) {
	permF <- sprintf("%s/markerMAF_perm.txt", outDir)
	perm_stats <- list.files(path=permMAF_F, pattern="[0-9].txt$", full.names=T)

	# Pool all permuted stats
	cat("Pooling perm .txt\n")
	system(sprintf("cat /dev/null > %s.tmp", permF))
	for (f in perm_stats) {
		print(f)
		system(sprintf("cat %s | cut -f 2 > %s.tmp", f, f)) #get second columns
	}																											#with deltaMAF values

	inDir <- dirname(perm_stats[1])
	system(sprintf("paste -d, %s/*.tmp > %s/perm_stats.txt", inDir, outDir))

	# Now write markerMAF_perm.txt
	system(sprintf("cat %s | cut -f1 > %s/marker.tmp", realF, outDir))
	system(sprintf("paste %s/marker.tmp %s/perm_stats.txt > %s.tmp",
								 outDir, outDir, permF))
	system(sprintf('sed "1s/.*/Marker\tCHI2_PERM/" %s.tmp > %s', permF, permF))

	# Remove all .tmp and redundant files
	unlink(sprintf("%s/*.tmp", permMAF_F))
	unlink(sprintf("%s/*.tmp", outDir))
	unlink(sprintf("%s/perm_stats.txt", outDir))

	# Run phenotype-permuted GSEA
	cat("* Running phenotype-based permutations.\n\n")
	cmd <- sprintf(phenoPerm, CALC_GSEA, realF, pathFile, snp2geneF,
								 permF, statOut, leOut, minGene, maxGene,
								 setSeed, snp2genedist, outDir)
	system(cmd)
} else {
	# Run genotype-permuted GSEA
	cat("* Running genotype-based permutations.\n\n")
	cmd <- sprintf(genoPerm, CALC_GSEA, realF, pathFile, snp2geneF,
								 statOut, leOut, minGene, maxGene,
								 setSeed, snp2genedist, outDir)
	system(cmd)
}

# Combine GSEA results
system(sprintf("perl %s %s/calculate_gsea.log > %s/combined_res.log",
							 COMB_GSEA, outDir, outDir))

# Format results for output
cat("Formatting results for output...")
dat <- read.delim(sprintf("%s/combined_res.log", outDir),skip=1,h=F,as.is=T)
dat[,1] <- sub("Geneset=","",	dat[,1])
dat[,2] <- sub("Size=","",		dat[,2])
dat[,3] <- sub("ES=","",		dat[,3])
dat[,4] <- sub("NES=","",		dat[,4])
dat[,5] <- sub("NominalP=","",	dat[,5])
dat[,6] <- sub("FDR=","",		dat[,6])
dat[,7] <- sub("FWER=","",		dat[,7])
colnames(dat) <- c("Geneset", "Size", "ES", "NES", "NominalP", "FDR", "FWER")
dat <- dat[order(dat$FDR),]
write.table(dat, file=resOut, sep="\t", col=T, row=F, quote=F)
cat(" done.\n")

# Write output to excel format
#write.xlsx(dat,file=sprintf("%s/results.xlsx", outDir), col=T, row=F)

# Tally number of pathways at sig. threshold FDR < 0.05 and FDR < 0.1
sigPaths1 <- which(dat$FDR < 0.1)
sigPaths2 <- which(dat$FDR < 0.05)
cat("\n------SIGNIFICANT GSEA RESULTS------\n")
cat("Number of significantly enriched pathways with FDR < 0.1 =",
		length(sigPaths1))
cat("\nNumber of significantly enriched pathways with FDR < 0.05 =",
		length(sigPaths2))

# Write table of significant pathways to output directory
write.table(dat[sigPaths1,], file=sprintf("%s/pathways_FDR0.1.txt",
						outDir), sep="\t", col=T, row=F, quote=F)
write.table(dat[sigPaths2,], file=sprintf("%s/pathways_FDR0.05.txt",
						outDir), sep="\t", col=T, row=F, quote=F)

# Format results for Cytoscape enrichment map
cat("\n\nFormatting output for use in Cytoscape...")
dat <- dat %>% mutate( Description = Geneset )
dat <- subset(dat, select=c(Geneset, NominalP, FDR, Description))
dat <- dat[c("Geneset", "Description", "NominalP", "FDR")]
write.table(dat, file=sprintf("%s/results_forEM.txt", outDir),
						sep="\t", col=T, row=F, quote=F)
cat(" done.\n")

# Cleanup
unlink(sprintf("%s/combined_res.log", outDir))
