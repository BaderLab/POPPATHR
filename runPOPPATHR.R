#!/usr/bin/env rscript

#-------------------------------------------------------------------------------
# POPPATHR: Population-based pathway analysis of SNP-SNP coevolution
# Special thanks to Shraddha Pai
# Last modified 26 August 2019

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
#system("./sh/annotations.sh")

geneF <- sprintf("%s/refGene.hg19.header.txt", annoDir) # NOTE originally downloaded 2017-10-30
pathF <- list.files(pattern="*.gmt", path=annoDir, full.names=TRUE)

#-------------------------------------------------------------------------------
## SOFTWARE DOWNLOAD ##
message("Downloading external software\n")
softwareDir <- sprintf("%s/software", dataDir)
if (!file.exists(softwareDir)) dir.create(softwareDir)

# Execute script to download external software
#system("./sh/software.sh")

#-------------------------------------------------------------------------------
## DATA DOWNLOAD ##
message("Downloading HapMap3 genotype data\n")
genoDir <- sprintf("%s/genotypes", dataDir)
if (!file.exists(genoDir)) dir.create(genoDir)

# Execute script to download SNP genotype data
#system("./sh/genotypes.sh")

genoF <- sprintf("%s/HM3_2010_05_phase3", genoDir)
popsF <- sprintf("%s/relationships_w_pops_041510.txt", genoDir)

#-------------------------------------------------------------------------------
## EXTERNAL R PACKAGES ##
library(librarian)
pkgs <- c("plyr", "dplyr", "ggplot2", "data.table", "stringr",
          "reshape2", "gdata", "RColorBrewer", "gridExtra", "cowplot",
          "GenomicRanges", "snpStats", "RCy3")
shelf(pkgs, cran_repo="https://cran.r-project.org")

# Load PopulationPathways R functions
Rfun <- list.files(pattern="*.R", path="R", full.names=TRUE)
mapply(source, Rfun)

#-------------------------------------------------------------------------------
## USER-DEFINED PARAMETERS ##
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
## STEP 1: Identify selection-enriched pathway gene sets
selPaths <- function(pop1, pop2) {

  ### Create dated output directory
  message("\n**Creating output directory.\n")
  dt <- format(Sys.Date(), "20%y%m%d")
  outDir <- sprintf("%s/%s_out_%s_%s", dataDir, dt, pop1, pop2)
  if (!file.exists(outDir)) dir.create(outDir)

  #sink(sprintf("%s/runPipeline.log", outDir))
  cat("Running pipeline with the following input data and parameters...\n")
  cat(sprintf("\tPopulations: %s and %s
               \tOutput directory: %s
               \tGenotyping dataset: %s
               \tPathway size limit: %s - %s genes
               \tSNP-gene mapping threshold: %skb
               \tPermutation cycles: %s\n",
      pop1, pop2, outDir, genoF, minGene, maxGene, snp2genedist/1000, setPerm))
  Sys.sleep(3)

  ### Recode PLINK fam file to population coding
  message("\n**Getting case/control (i.e., population) status.\n")
  # New PLINK file names after population recoding
  famName <- sprintf("%s_%s_%s", basename(genoF), pop1, pop2)
  recodeFAM(genoF=genoF, pop1=pop1, pop2=pop2, popsF=popsF, outF=famName)
  Sys.sleep(3)

  ### Determine population substructure via PLINK --cluster
  message(sprintf("\n**Determining population substructure.\n"))
  pcaDir <- sprintf("%s/pca", outDir)
  if (!file.exists(pcaDir)) dir.create(pcaDir)

  realFam <- sprintf("%s_%s_%s.fam", genoF, pop1, pop2)
  pcaF <- sprintf("%s/%s", pcaDir, famName)
  popPCA(genoF=genoF, famF=realFam, pop1=pop1, pop2=pop2, outF=pcaF)
  Sys.sleep(3)

  ### Calculate FST estimation per SNP between both populations
  message("\n**Calculating population SNP-level FST.\n")
  fstDir <- sprintf("%s/fst", outDir)
  if (!file.exists(fstDir)) dir.create(fstDir)

  fstF <- sprintf("%s/markerFST.txt", fstDir)
  calcFST(genoF=genoF, realFam=realFam, outF=fstF, outDir=fstDir)
  Sys.sleep(3)

  ### Map input SNPs to genes
  message("\n**Mapping input SNPs to genes.\n")
  gseaDir <- sprintf("%s/gsea", outDir)
  if (!file.exists(gseaDir)) dir.create(gseaDir)

  snpF <- sprintf("%s.bim", genoF)
  snp2geneF <- sprintf("%s/snp2gene.txt", gseaDir)
  SNP2gene(inF=snpF, geneF=geneF, outF=snp2geneF)
  Sys.sleep(3)

  ### Run GSEA
  message("\n**Running gene-set enrichment analysis.\n")
  setupGSEArun(realF=fstF, pathF=pathF,
               snp2geneF=snp2geneF, snp2genedist=snp2genedist,
               minGene=minGene, maxGene=maxGene,
               setPerm=setPerm, outDir=gseaDir)
  Sys.sleep(3)
}

# Run CEU vs. YRI and CEU vs. LWK analyses
pops_1 <- c("CEU", "CEU")
pops_2 <- c("YRI", "LWK")
mapply(selPaths, pops_1, pops_2)

#-------------------------------------------------------------------------------
## STEP 2: Identify genetic coevolution within and bewteen pathways
### Generate SNP lists per selection-enriched and unenriched pathway
message("\n**Generating SNP lists per selection-enriched and unenriched pathway(s).\n")

### Write into directory of first analysis (eg. CEU_YRI)
pop1=pops_1[1]
pop2=pops_2[1]

outDir <- list.files(pattern=sprintf("*%s_%s", pop1, pop2), path=getwd())
ldDir       <- sprintf("%s/ld", outDir)
enrichDir   <- sprintf("%s/enriched_paths", ldDir)
unenrichDir <- sprintf("%s/unenriched_paths", ldDir)

if (!file.exists(ldDir)) dir.create(ldDir)
if (!file.exists(enrichDir)) dir.create(enrichDir)
if (!file.exists(unenrichDir)) dir.create(unenrichDir)

realFam <- sprintf("%s_%s_%s.fam", genoF, pop1, pop2)

### Get GSEA results files from each analysis
resDirs <- list.files(pattern="*out", path=getwd())
resF <- list.files(pattern="results.txt", path=sprintf("%s/gsea", resDirs), full.names=TRUE)
gseaStatF <- sprintf("%s/gsea/gseaStatFile.txt", outDir)
snp2geneF <- sprintf("%s/gsea/snp2gene.txt", outDir)

getPathStats(genoF=genoF, resF=resF,
             gseaStatF=gseaStatF, snp2geneF=snp2geneF,
             enrichNES=enrichNEScut, unenrichNES=unenrichNEScut,
						 enrichDir=enrichDir, unenrichDir=unenrichDir)

### Calculate LD statisics for each enriched pathway compared
message("\n**Calculating trans-chromosomal LD statistics.\n")

#WPM (within-pathway model)
wpmDir <- sprintf("%s/WPM", ldDir)
if (!file.exists(wpmDir)) dir.create(wpmDir)
LDstatsWPM(enrichDir=enrichDir, unenrichDir=unenrichDir,
           pop1=pop1, pop2=pop2, outDir=wpmDir)

#BPM (between-pathway model)
bpmDir <- sprintf("%s/BPM", ldDir)
if (!file.exists(bpmDir)) dir.create(bpmDir)
LDstatsBPM(enrichDir=enrichDir, unenrichDir=unenrichDir,
           pop1=pop1, pop2=pop2, snp2geneF=snp2geneF, outDir=bpmDir)

#-------------------------------------------------------------------------------
## STEP3: Getting gene properties for selection-enriched genes / variants
message("\n**Getting selection-enriched gene properties.\n")

propDir <- sprintf("%s/geneProp", outDir)
if (!file.exists(propDir)) dir.create(propDir)

geneDir <- enrichDir
geneF  <- sprintf("%s/genes_leadingEdge_unique_hc.txt", geneDir)
clustF <- sprintf("%s/gseaStat_per_hc_pathway_updated.txt", geneDir)

outDir=propDir
