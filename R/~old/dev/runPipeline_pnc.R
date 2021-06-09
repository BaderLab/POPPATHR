#' Runs pathway analysis with PLINK + GSEA using the difference in SNP
#' level FST between two populations (previously dMAF)

#-------------------------------------------------------------------------------
## INPUT PARAMETERS
dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways"
#dataDir <- "/Users/catherineross/PopulationPathways"

data <- "PNC"
pop1 <- "CEU"
pop2 <- "ASW"
geno <- "PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW"

# Run GSEA using phenotype permutations - yes (TRUE) or no (FALSE)
# If FALSE, will use genotype permutations (see setupGSEArun.R for details)
phenoPerm <- TRUE

  # If phenoPerm = TRUE, set number of permutations to run + random seed
  numPerm   <- 1000L
  setSeed   <- 42L

# Indicate which SNP-level test statistic to compute, MAF or FST (03/24/2018)
snpStat <- "FST" #FST or MAF

# Use absolute ΔMAF as SNP statistic - yes (TRUE) or no (FALSE; will use ΔMAF)?
absMAF <- TRUE

# Set SNP mapping args
filterSNP  <- FALSE

  # If filterSNP = TRUE, set max. distance away from gene to filter out SNPs
  # (i.e., if filterDist = 20000, filters out all SNPs more than 20,000 bp
  # away from nearest gene)
  filterDist <- 20000L

# Set GSEA args
setPerm      <- 10000L
minGene      <- 10L
maxGene      <- 300L
snp2genedist <- 500000L # testing wider SNP-gene mapping filter
#minGene       <- 20L
#maxGene       <- 200L
#snp2genedist  <- 10000L

#-------------------------------------------------------------------------------
## SOFTWARE & ANNOTATION
PLINK 		<- "/media/catherine/DATAPART1/Software/plink_linux_1.90/plink"
CALC_GSEA <- "/media/catherine/DATAPART1/Software/GenGen-1.0.1/calculate_gsea.pl"
COMB_GSEA	<- "/media/catherine/DATAPART1/Software/GenGen-1.0.1/combine_gsea.pl"
#PLINK <- "/Users/catherineross/Software/plink_mac/plink"
#CALC_GSEA <- "/Users/catherineross/Software/GenGen-1.0.1/calculate_gsea.pl"
#COMB_GSEA <- "/Users/catherineross/Software/GenGen-1.0.1/combine_gsea.pl"

geneFile <- sprintf("%s/anno/refGene/refGene.hg19.header.txt", dataDir)
pathFile <- sprintf(paste0("%s/anno/baderlab/Human_GOBP_AllPathways_no_GO",
                           "_iea_April_24_2016_symbol.gmt"), dataDir)
#pathFile <- sprintf(paste0("%s/anno/baderlab/Human_GOBP_AllPathways_no_GO_iea_",
#                           "April_24_2016_symbol_sans_HLA_IFN.gmt"), dataDir)
#pathFile <- sprintf(paste0("%s/anno/baderlab/Human_GOBP_AllPathways_no_GO",
#                           "_iea_April_24_2016_symbol_selection_genes.gmt"), dataDir)

## EXTERNAL PACKAGES
suppressMessages({
  require(plyr)            #LDstats
  require(dplyr)           #calcMAFdiff, getSNPlists, LDstats, pathPosSel
  require(ggplot2)         #popPCA, LDstats, pathPosSel
  require(GenomicRanges)   #SNP2gene
  require(data.table)      #SNP2gene, calcMAFdiff
  require(tools)           #calcMAFdiff, getPathStats
  require(stringr)         #getPathStats
  require(snpStats)        #LDstats
  require(reshape2)        #LDstats, pathPosSel
  require(RColorBrewer)    #LDstats
  require(gridExtra)       #LDstats
  require(cowplot)         #LDstats
  require(splitstackshape) #pathPosSel
  require(utils)           #interactive
  require(SNPRelate)       #LDstats (composite LD)
})

#-------------------------------------------------------------------------------
## WORK BEGINS
dt	    <- format(Sys.Date(), "%y%m%d")
#outDir  <- sprintf("%s/res/out_%s_%s_FST_%s-%s_20-200gene",
#                   dataDir, dt, geno, pop1, pop2)
outDir  <- sprintf("%s/res/out_%s_%s_%s_%s-%s_%s-%sgene_%skb-dist_1000pheno-perm",
                   dataDir, dt, geno, snpStat, pop1, pop2, minGene, maxGene,
                   snp2genedist/1000)

if (file.exists(outDir)) {
  cat("Directory exists! Not overwriting\n")
  Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)

#-------------------------------------------------------------------------------
sink(sprintf("%s/runPipeline.log", outDir))
cat("Running pipeline with the following input data and parameters...\n")
cat(sprintf("\tPopulations: %s and %s
             \tGenotyping dataset: %s
             \tSNP-level test statistic: %s
             \tPathway size limit: %s - %s genes
             \tSNP-gene mapping threshold: %skb
             \tPermutation cycles: %s\n",
      pop1, pop2, data, snpStat, minGene, maxGene, snp2genedist/1000, numPerm))

message("\nSTEP 1: Getting case/control (i.e., population) status.\n")

genoF   <- sprintf("%s/data/%s/all/%s", dataDir, data, geno)
realFAM <- sprintf("%s.fam", genoF)

source("recodeFAM.R")
recodeFAM(genoF, popInfo, recode=FALSE,
          phenoPerm=phenoPerm, numPerm=numPerm, setSeed)

## Calculate difference in minor allele frequency (MAF) for SNPs
## between both populations
message(sprintf("\nCalculating SNP test statistic (Δ%s).\n", snpStat))
plinkDir <- sprintf("%s/freq", outDir)
if (!file.exists(plinkDir)) dir.create(plinkDir)

if (snpStat == "MAF") {
  source("calcMAFdiff.R")
  calcMAFdiff(genoF, realFAM, phenoPerm, permFAM=permDir,
              PLINK=PLINK, absMAF=absMAF, outDir=plinkDir)
} else if (snpStat == "FST") {
  source("calcFST.R")
  calcFST(genoF, realFAM, PLINK=PLINK,
          phenoPerm, permFAM=sprintf("%s/perm", outDir),
          outDir=plinkDir)
}
#-------------------------------------------------------------------------------
## Map input SNPs to genes
message("\nMapping input SNPs to genes.\n")
gseaDir <- sprintf("%s/gsea", outDir)
if (!file.exists(gseaDir)) dir.create(gseaDir)

source("SNP2gene.R")
SNP2gene(snpFile=sprintf("%s.bim", genoF), geneFile,
         filterSNP=filterSNP,
         outFile=sprintf("%s/snp2gene.txt", gseaDir))

#-------------------------------------------------------------------------------
## Run GSEA
message("\nRunning gene-set enrichment analysis.\n")
resDir  <- sprintf("%s/pathway_analysis", gseaDir)
if (!file.exists(resDir)) dir.create(resDir)

snp2geneF <- sprintf("%s/snp2gene.txt", gseaDir)
realF   <- sprintf("%s/marker%s.txt", plinkDir, snpStat)

source("setupGSEArun.R")
setupGSEArun(realF, phenoPerm=phenoPerm, permFST_F=sprintf("%s/perm", outDir),
             pathFile, snp2geneF, snp2genedist=snp2genedist,
             CALC_GSEA=CALC_GSEA, COMB_GSEA=COMB_GSEA,
             minGene=minGene, maxGene=maxGene,
             setPerm=setPerm, setSeed=setSeed, outDir=resDir)

sink()
