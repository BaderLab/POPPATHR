#' Runs pathway analysis with PLINK + GSEA using the difference in SNP
#' level FST between two populations

#-------------------------------------------------------------------------------
## INPUT PARAMETERS
dataDir <- "/Users/catherineross/PopulationPathways"      # parent directory
genoDir <- sprintf("%s/data/HM3_2010-05_phase3", dataDir) # genotype input directory

# Specify two populations to test
pop1 <- "CEU"
pop2 <- "YRI"

# Get genotype + population information
geno <- list.files("*.bed", path=genoDir, full.names=T)
if (!length(geno)) {
  stop("PLINK-formatted genotype file not found! (ie. bed, bim, fam files)")
} else {
  genoF <- substr(basename(geno), 0, nchar(basename(geno))-4)
  popsF <- list.files("relationships*", path=genoDir, full.names=T)
}

# Run GSEA using phenotype permutations - yes (TRUE) or no (FALSE)
# If FALSE, will use genotype permutations (see setupGSEArun.R for details)
phenoPerm <- FALSE

  # If phenoPerm = TRUE, set number of permutations to run + random seed
  numPerm   <- 100L
  setSeed   <- 42L

# Indicate which SNP-level test statistic to compute
snpStat <- "FST" # previously supported using MAF

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
snp2genedist <- 500000L  # testing wider SNP-gene mapping filter

# GSEA enrichment filtering args
LEsubset <- FALSE
lcNEScut <- 0.1   # NES threshold to define unenriched pathways
hcFDRcut <- 0.05  # FDR threshold to define selection-enriched pathways

#-------------------------------------------------------------------------------
## SOFTWARE & ANNOTATION
PLINK <- "/Users/catherineross/Software/plink_mac/plink"
CALC_GSEA <- "/Users/catherineross/Software/GenGen-1.0.1/calculate_gsea.pl"
COMB_GSEA <- "/Users/catherineross/Software/GenGen-1.0.1/combine_gsea.pl"

geneFile <- sprintf("%s/anno/refGene/refGene.hg19.header.txt", dataDir)
pathFile <- sprintf(paste0("%s/anno/baderlab/Human_GOBP_AllPathways_no_GO",
                           "_iea_April_24_2016_symbol.gmt"), dataDir)

## EXTERNAL PACKAGES
## NOTE use librarian to install / attach
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
  require(gdata)           #LDstats
  require(RColorBrewer)    #LDstats
  require(gridExtra)       #LDstats
  require(cowplot)         #LDstats
  require(utils)           #interactive
})

# Load publication-style ggplot parameters
#source("themePublication.R")
#-------------------------------------------------------------------------------
## WORK BEGINS
dt	    <- format(Sys.Date(), "%y%m%d")
outDir  <- sprintf("%s/res/out_%s_%s_%s-%s_%s_%s-%sgene_%skb-dist_%sperm_updated_data",
                   dataDir, dt, genoF, pop1, pop2, snpStat, minGene, maxGene,
                   snp2genedist/1000, setPerm)

if (file.exists(outDir)) {
  cat("Directory exists! Not overwriting\n")
  Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)

#-------------------------------------------------------------------------------
## STEP 1: Recode PLINK .fam file to case/control coding
sink(sprintf("%s/runPipeline.log", outDir))
cat("Running pipeline with the following input data and parameters...\n")
cat(sprintf("\tPopulations: %s and %s
             \tGenotyping dataset: %s
             \tSNP-level test statistic: %s
             \tPathway size limit: %s - %s genes
             \tSNP-gene mapping threshold: %skb
             \tPermutation cycles: %s\n",
      pop1, pop2, genoF, snpStat, minGene, maxGene, snp2genedist/1000, setPerm))

message("\nSTEP 1: Getting case/control (i.e., population) status.\n")

source("recodeFAM.R")
recodeFAM(genoF, pop1, pop2, popsF, recode=TRUE,
          phenoPerm=phenoPerm, numPerm=numPerm, setSeed,
          famName=sprintf("%s_%s_%s", geno, pop1, pop2))

#-------------------------------------------------------------------------------
## STEP 2: Determine population substructure via PLINK --cluster
message(sprintf("\nSTEP 2: Determining population substructure.\n"))
pcaDir <- sprintf("%s/pca", outDir)
if (!file.exists(pcaDir)) dir.create(pcaDir)

## new recoded fam file (case/control coded)
realFAM <- sprintf("%s_%s_%s.fam", genoF, pop1, pop2)
ggtitle <- "Population stratification via PCA"

source("popPCA.R")
#popPCA(genoF, famF=realFAM, dimensions=3L, PLINK=PLINK,
#       outF=sprintf("%s/%s_%s", pcaDir, pop1, pop2))

#-------------------------------------------------------------------------------
## STEP 3: Calculate FST estimation per SNP
## between both populations
message(sprintf("\nSTEP 3: Calculating SNP test statistic (%s).\n", snpStat))
plinkDir <- sprintf("%s/freq", outDir)
if (!file.exists(plinkDir)) dir.create(plinkDir)

source("calcFST.R")
calcFST(genoF, realFAM, PLINK=PLINK,
        phenoPerm, permFAM=sprintf("%s/perm", outDir),
        outDir=plinkDir)

#-------------------------------------------------------------------------------
## STEP 4: Map input SNPs to genes
message("\nSTEP 4: Mapping input SNPs to genes.\n")
gseaDir <- sprintf("%s/gsea", outDir)
if (!file.exists(gseaDir)) dir.create(gseaDir)

source("SNP2gene.R")
SNP2gene(snpFile=sprintf("%s.bim", genoF), geneFile,
         filterSNP=filterSNP,
         outFile=sprintf("%s/snp2gene.txt", gseaDir))

#-------------------------------------------------------------------------------
## STEP 5: Run GSEA
message("\nSTEP 5: Running gene-set enrichment analysis.\n")
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

#-------------------------------------------------------------------------------
## STEP 6: Generate SNP lists per selection-enriched and unenriched pathway
## NOTE: DO THIS STEP AFTER RUNNING TWO POP COMPARISONS ##
message("\nSTEP 6: Generating SNP lists per selection-enriched and unenriched pathway(s).\n")
ldDir   <- sprintf("%s/ld", outDir)
hcSnps  <- sprintf("%s/hc_snps", ldDir)   # selection-enriched
lcSnps  <- sprintf("%s/lc_snps", ldDir) # unenriched

if (!file.exists(ldDir)) dir.create(ldDir)
if (!file.exists(hcSnps)) dir.create(hcSnps)
if (!file.exists(lcSnps)) dir.create(lcSnps)

## NOTE: interactive way to select replication directory
switch(menu(c("Yes", "No"), title="Have you run a replication analysis?"),
    cat("Next.\n"),
    cat("Exiting pipeline.\n"))

resDirs <- c(outDir, repDir)
resF <- sprintf("%s/%s/%s/results.txt", resDirs, basename(gseaDir), basename(resDir))

gseaStatF <- sprintf("%s/gseaStatFile.txt", dirname(resF))
gseaLEoutF <- sprintf("%s/gseaLEout.txt", dirname(resF))
testStatF <- sprintf("%s/%s/%s/marker%s.txt",
                     resDirs, basename(gseaDir), basename(plinkDir), snpStat)

source("getPathStats.R")
getPathStats(genoF, resF, gseaStatF, snp2geneF, PLINK,
             LEsubset=FALSE, lcNEScut=lcNEScut, hcFDRcut=hcFDRcut,
						 hcOutDir=hcSnps, lcOutDir=lcSnps)

#-------------------------------------------------------------------------------
## STEP 7: Calculate LD statisics for each enriched pathway compared
## to the cumulative set of nonenriched pathways
message("\nSTEP 8: Calculating inter-chromosomal LD statistics.\n")
# Define statistic for measuring association between SNP pairs
statistic <- "R.squared" # or D.prime

#WPM (within-pathway model)
statDir <- sprintf("%s/res_%s_%s_WPM", ldDir, basename(hcSnps), basename(lcSnps))
if (!file.exists(statDir)) dir.create(statDir)

source("LDstatsWPM.R")
LDstats(hcInDir=hcSnps, lcInDir=lcSnps, statistic=statistic,
        popNames=popNames, outDir=statDir)

#BPM (between-pathway model)
statDir <- sprintf("%s/res_%s_%s_BPM", ldDir, basename(hcSnps), basename(lcSnps))
if (!file.exists(statDir)) dir.create(statDir)

source("LDstatsBPM.R")
LDstats(hcInDir=hcSnps, lcInDir=lcSnps, popNames=popNames,
        snp2geneF=snp2geneF, outDir=statDir)
