#!/usr/bin/env rscript

#-------------------------------------------------------------------------------
# POPPATHR: Population-based pathway analysis of SNP-SNP coevolution
# Special thanks to Shraddha Pai
# Last modified 8 August 2019

#-------------------------------------------------------------------------------
## PREPARE PIPELINE ##
message("Preparing 'PopPaths' pipeline!")
system("chmod +x sh/*.sh")
dataDir <- getwd()

#-------------------------------------------------------------------------------
## ANNOTATION FILES ##
message("Downloading genome and pathway annotation files\n")
annoDir <- sprintf("%s/anno", dataDir)
if (!file.exists(annoDir)) dir.create(annoDir)

# Execute script to get annotation files
system("./sh/annotations.sh")

geneF <- sprintf("%s/refGene.hg19.header.txt", annoDir) # NOTE originally downloaded 2017-10-30
pathF <- list.files(pattern="*.gmt", path=annoDir, full.names=TRUE)

#-------------------------------------------------------------------------------
## SOFTWARE DOWNLOAD ##
message("Downloading external software\n")
softwareDir <- sprintf("%s/software", dataDir)
if (!file.exists(softwareDir)) dir.create(softwareDir)

# Execute script to download external software
system("./sh/software.sh")

#-------------------------------------------------------------------------------
## DATA DOWNLOAD ##
message("Downloading HapMap3 genotype data\n")
genoDir <- sprintf("%s/genotypes", dataDir)
if (!file.exists(genoDir)) dir.create(genoDir)

# Execute script to download SNP genotype data
system("./sh/genotypes.sh")

genoF <- sprintf("%s/HM3_2010_05_phase3", genoDir)
popsF <- sprintf("%s/relationships_w_pops_041510.txt", genoDir)

#-------------------------------------------------------------------------------
## EXTERNAL R PACKAGES ##
library(librarian)
pkgs <- c("plyr", "dplyr", "ggplot2", "data.table", "tools", "stringr",
          "reshape2", "gdata", "RColorBrewer", "gridExtra", "cowplot",
          "utils", "GenomicRanges", "snpStats")
shelf(pkgs, cran_repo="https://cran.r-project.org")

# Load PopulationPathways R functions
Rfun <- list.files(pattern="*.R", path="R", full.names=TRUE)
mapply(source, Rfun)

#-------------------------------------------------------------------------------
## USER-DEFINED PARAMETERS ##
# Tested populations
pop1 <- "CEU"
pop2 <- "YRI"

# GSEA function args
setPerm      <- 10000L  # Number of permutation cycles
minGene      <- 10L     # Min. number of genes in pathway
maxGene      <- 300L    # Max. number of genes in pathway
snp2genedist <- 500000L # Max. SNP-gene mapping distance

# GSEA result filtering args
enrichNES   <- 3   # NES threshold to define selection-enriched pathways
unenrichNES <- 0.1 # abs NES threshold to define unenriched pathways

#-------------------------------------------------------------------------------
## WORK BEGINS ##
dt <- format(Sys.Date(), "20%y%m%d")
outDir <- sprintf("%s/%s_out", dataDir, dt)
if (!file.exists(outDir)) dir.create(outDir)

sink(sprintf("%s/runPipeline.log", outDir))
cat("Running pipeline with the following input data and parameters...\n")
cat(sprintf("\tPopulations: %s and %s
             \tGenotyping dataset: %s
             \tPathway size limit: %s - %s genes
             \tSNP-gene mapping threshold: %skb
             \tPermutation cycles: %s\n",
      pop1, pop2, genoF, minGene, maxGene, snp2genedist/1000, setPerm))

#-------------------------------------------------------------------------------
## STEP 1: Recode PLINK fam file to population coding
message("\nSTEP 1: Getting case/control (i.e., population) status.\n")

# New PLINK file names after population recoding
famName <- sprintf("%s_%s_%s", basename(genoF), pop1, pop2)
recodeFAM(genoF=genoF, pop1=pop1, pop2=pop2, popsF=popsF, outF=famName)

#-------------------------------------------------------------------------------
## STEP 2: Determine population substructure via PLINK --cluster
message(sprintf("\nSTEP 2: Determining population substructure.\n"))
pcaDir <- sprintf("%s/pca", outDir)
if (!file.exists(pcaDir)) dir.create(pcaDir)

realFam <- sprintf("%s_%s_%s.fam", genoF, pop1, pop2)
pcaF <- sprintf("%s/%s", pcaDir, famName)
popPCA(genoF=genoF, famF=realFam, outF=pcaF)

#-------------------------------------------------------------------------------
## STEP 3: Calculate FST estimation per SNP
## between both populations
message("\nSTEP 3: Calculating population SNP-level FST.\n")
fstDir <- sprintf("%s/fst", outDir)
if (!file.exists(fstDir)) dir.create(fstDir)

fstF <- sprintf("%s/markerFST.txt", fstDir)
calcFST(genoF=genoF, realFam=realFam, outF=fstF, outDir=fstDir)

#-------------------------------------------------------------------------------
## STEP 4: Map input SNPs to genes
message("\nSTEP 4: Mapping input SNPs to genes.\n")
gseaDir <- sprintf("%s/gsea", outDir)
if (!file.exists(gseaDir)) dir.create(gseaDir)

snpF <- sprintf("%s.bim", genoF)
snp2geneF <- sprintf("%s/snp2gene.txt", gseaDir)
SNP2gene(inF=snpF, geneF=geneF, outF=snp2geneF)

#-------------------------------------------------------------------------------
## STEP 5: Run GSEA
message("\nSTEP 5: Running gene-set enrichment analysis.\n")
setupGSEArun(realF=fstF, pathF=pathF,
             snp2geneF=snp2geneF, snp2genedist=snp2genedist,
             minGene=minGene, maxGene=maxGene,
             setPerm=setPerm, outDir=gseaDir)



#-------------------------------------------------------------------------------
## STEP 6: Generate SNP lists per selection-enriched and unenriched pathway
message("\nSTEP 6: Generating SNP lists per selection-enriched and unenriched pathway(s).\n")
ldDir <- sprintf("%s/ld", outDir)
enrich   <- sprintf("%s/enriched", ldDir)
unenrich <- sprintf("%s/unenriched", ldDir)

if (!file.exists(ldDir)) dir.create(ldDir)
if (!file.exists(enrich)) dir.create(enrich)
if (!file.exists(unenrich)) dir.create(unenrich)

resDirs <- c(outDir, repDir)
resF <- sprintf("%s/%s/%s/results.txt", resDirs, basename(gseaDir), basename(resDir))

gseaStatF <- sprintf("%s/gseaStatFile.txt", dirname(resF))
gseaLEoutF <- sprintf("%s/gseaLEout.txt", dirname(resF))
testStatF <- sprintf("%s/%s/%s/marker%s.txt",
                     resDirs, basename(gseaDir), basename(plinkDir), snpStat)

source("getPathStats.R")
getPathStats(genoF=genoF, resF, gseaStatF, snp2geneF,
             enrichNES=hcNEScut, unenrichNES=lcNEScut,
						 enrichDir=enrich, unenrichDir=unenrich)

#-------------------------------------------------------------------------------
## STEP 7: Calculate LD statisics for each enriched pathway compared
## to cumulative set of unenriched pathways
message("\nSTEP 7: Calculating inter-chromosomal LD statistics.\n")
# Define statistic for measuring association between SNP pairs
statistic <- "R.squared"

#WPM (within-pathway model)
statDir <- sprintf("%s/%s_%s_WPM", ldDir, basename(enrich), basename(unenrich))
if (!file.exists(statDir)) dir.create(statDir)

source("LDstatsWPM.R")
LDstats(hcInDir=enrich, lcInDir=unenrich, statistic=statistic,
        popNames=popNames, outDir=statDir)

#BPM (between-pathway model)
statDir <- sprintf("%s/%s_%s_BPM", ldDir, basename(enrich), basename(unenrich))
if (!file.exists(statDir)) dir.create(statDir)

source("LDstatsBPM.R")
LDstats(hcInDir=enrich, lcInDir=unenrich, popNames=popNames,
        snp2geneF=snp2geneF, outDir=statDir)

#-------------------------------------------------------------------------------
## STEP 8: Getting gene properties for selection-enriched genes / variants
message("\nSTEP 8: Getting selection-enriched gene properties.\n")

propDir <- sprintf("%s/geneProp", outDir)
if (!file.exists(propDir)) dir.create(propDir)

geneDir <- hcSnps
geneF  <- sprintf("%s/genes_leadingEdge_unique_hc.txt", geneDir)
clustF <- sprintf("%s/gseaStat_per_hc_pathway_updated.txt", geneDir)

outDir=propDir
