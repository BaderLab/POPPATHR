#' runs pathway analysis with plink + GSEA

#-------------------------------------------------------------
# Input files

### CROSS replace with root to your WTCCC folder
dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways"
geno<- sprintf("%s/data/WTCCC/CD_case_control.autosomes.FINAL.grch38.sorted",dataDir)

dt		<- format(Sys.Date(),"%y%m%d")
outDir	<- sprintf("%s/results/out_%s",dataDir,dt)

#-------------------------------------------------------------
# Software and annotation
PLINK 		<- "/media/catherine/DATAPART1/Software/plink_linux_1.90/plink"
CALC_GSEA <- "/media/catherine/DATAPART1/Software/GenGen-1.0.1/calculate_gsea.pl"
COMB_GSEA	<- "/media/catherine/DATAPART1/Software/GenGen-1.0.1/combine_gsea.pl"

geneFile <- sprintf("%s/anno/refGene/refGene.hg19.header.txt", dataDir)
pathFile <- sprintf(paste0("%s/anno/baderlab/Human_GOBP_AllPathways_no_GO",
                           "_iea_April_24_2016_symbol.gmt"), dataDir)

#--------------------------------------------------------------
# Work begins
if (file.exists(outDir)) {
	cat("Directory exists! Not overwriting\n")
	Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)

require(GWAS2pathway)
plinkDir <- sprintf("%s/gwas",outDir)
if (!file.exists(plinkDir)) dir.create(plinkDir)

# --------------------------------------
# GWAS

# create the plink jobs for the real GWAS and 100 GWAS where
# case control labels have been permuted
jobInfo <- setupPlink(geno,numPerm=100L,outDir=plinkDir,PLINK=PLINK,
		numCores=2L,
		   jobFile=sprintf("%s/plinkJobs.txt",plinkDir))

# now run the plink jobs.
# note: you should have installed GNU parallel to use this
# on os x, install GNU parallel via homebrew.
# First get homebrew if you don't have it (http://brew.sh/)
# Then on bash:
# $ (sudo) brew install parallel

# examine the signal
#dat <- read.table(jobInfo$outFile,h=T,as.is=T)
#dat <- na.omit(dat)
#hist(dat$P)

## uncomment to look at Manhattan plot
#require(qqman)
#gwasRes <- dat[,c("SNP","CHR","BP","P")]
#manhattan(gwasRes,ylim=c(0,20),cex=0.4,col=c("purple","blue"))

# --------------------------------------
# Pathway analysis by GSEA

# map SNPs to genes
GSEAdir <- sprintf("%s/gsea", outDir)
if (!file.exists(GSEAdir)) dir.create(GSEAdir)
snpFile <- sprintf("%s/snp2gene.txt",GSEAdir)
map_SNP2gene(jobInfo$outFile,geneFile,snpFile)

# setup files for gsea pathway analysis
permFiles <- dir(path=sprintf("%s/perm",plinkDir),pattern=".assoc$")
permFiles <- paste(sprintf("%s/perm",plinkDir),permFiles,sep="/")
gseaFiles <- GSEA_setup(jobInfo$outFile,permFiles,GSEAdir)

# run GSEA
runGSEA(gseaFiles[1],gseaFiles[2],pathFile,snpFile,
		CALC_GSEA=CALC_GSEA, COMB_GSEA=COMB_GSEA)
