#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# POPPATHR: Population-based pathway analysis of SNP-SNP coevolution
# Special thanks to Shraddha Pai
# Last modified 28 February 2020
#------------------------------------------------------------------------------#

# Load all packages in a way that allows exporting to child environments
packages <- c("tidyverse", "data.table", "reshape2", "gdata", "RColorBrewer",
              "GenomicRanges", "fgsea", "argparse")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Load custom functions
Rfun <- list.files(pattern="*.R", path="R", full.names=TRUE)
for (file in Rfun) { source(file) }

######
# PARSE USER ARGUMENTS
######

# Make argparse object and arguments
parser <- ArgumentParser(
    description=paste("POPPATHR: Population-based pathway analysis of SNP-SNP coevolution.",
    "This script identifies selection-enriched pathways between two population cohorts.")
)
parser$add_argument("-v", "--verbose", action="store_false",
                    help="Print verbose output")
parser$add_argument("-p", "--populationPair", type="character",
                    help="Names of two population cohorts to test (e.g., CEU_YRI or CEU_LWK)")
parser$add_argument("-g", "--genotypeFile", type="character", default="data/genotypes/HM3_2010_05_phase3",
                    help="Path to PLINK (bed, bim, fam formatted) SNP genotype files [default %(default)s]")
parser$add_argument("-t", "--populationTable", type="character", default="data/genotypes/relationships_w_pops_041510.txt",
                    help="Path to table defining population genotypes [default %(default)s]")
parser$add_argument("-a", "--annotationFile", type="character", default="data/annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                    help="Path to gmt file containing gene set annotations [default %(default)s]")
parser$add_argument("-r", "--refgeneFile", type="character", default="data/annotations/refGene.hg19.header.txt",
                    help="Path to refGene genome annotation file [default %(default)s]")
parser$add_argument("-o", "--outputFolder", type="character", default="output",
                    help="Path to output folder [default %(default)s]")
parser$add_argument("--setPerm", type="integer", default=10000,
                    help="Number of GSEA permutation cycles to run [default %(default)s]")
parser$add_argument("--minGene", type="integer", default=10,
                    help="Minimum number of genes permitted in gene set [default %(default)s]")
parser$add_argument("--maxGene", type="integer", default=500,
                    help="Maximum number of genes permitted in gene set [default %(default)s]")
parser$add_argument("--scoreType", type="character", default="pos",
                    help="Defines GSEA score type [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

## NOTE LINES FOR TESTING
popPair <- "CEU_YRI"
genotypeFile <- "data/genotypes/HM3_2010_05_phase3"
populationTable <- "data/genotypes/relationships_w_pops_041510.txt"
annotationFile <- "data/annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt"
refgeneFile <- "data/annotations/refGene.hg19.header.txt"
outputFolder <- "output"
setPerm = 10000
minGene = 10
maxGene = 500
scoreType = "pos"
##

# Sets important parameters
popPair <- args$populationPair
genotypeFile <- args$genotypeFile
populationTable <- args$populationTable
annotationFile <- args$annotationFile
refgeneFile <- args$refgeneFile
outputFolder <- args$outputFolder
setPerm <- args$setPerm
minGene <- args$minGene
maxGene <- args$maxGene
scoreType <- args$scoreType

# Define I/O file names based on user input information
popOne <- unlist(lapply(strsplit(popPair, "_"), "[", 1))
popTwo <- unlist(lapply(strsplit(popPair, "_"), "[", 2))
snpFile <- sprintf("%s.bim", genotypeFile)
famName <- sprintf("%s_%s", basename(genotypeFile), popPair)
famFile <- sprintf("%s_%s.fam", genotypeFile, popPair)

# Checks if files exist
if (!file.exists(snpFile)) {
  stop("PLINK-formatted SNP genotype file does not exist!")
}
if (is.na(popOne)) {
  stop("Is --populationPair (-p) argument in '[POPULATION_NAME]_[POPULATION_NAME]' format?")
}
if (is.na(popTwo)) {
  stop("Is --populationPair (-p) argument in '[POPULATION_NAME]_[POPULATION_NAME]' format?")
}

# Create parent output folder if it doesn't already exist
if (!dir.exists(outputFolder)) { dir.create(outputFolder) }

# Create analysis-specific subfolders
popFolder <- sprintf("%s/%s", outputFolder, popPair)
if (!dir.exists(popFolder)) { dir.create(popFolder) }
pcaFolder <- sprintf("%s/pca", popFolder)
if (!file.exists(pcaFolder)) dir.create(pcaFolder)
fstFolder <- sprintf("%s/fst", popFolder)
if (!file.exists(fstFolder)) dir.create(fstFolder)
gseaFolder <- sprintf("%s/gsea", popFolder)
if (!file.exists(gseaFolder)) dir.create(gseaFolder)

# Pre-define downstream files
pcaFile <- sprintf("%s/%s", pcaFolder, famName)
fstFile <- sprintf("%s/markerFST.txt", fstFolder)
snp2geneFile <- sprintf("%s/snp2gene.txt", gseaFolder)

#------------------------------------------------------------------------------#
# IDENTIFY SELECTION-ENRICHED GENE SETS
#------------------------------------------------------------------------------#

######
# WORK BEGINS
######

# Sink cats to file
sink(sprintf("%s/getEnrichment.log", popFolder))

# Print user parameters to log
cat("POPPATHR population-based gene set enrichment analysis\n")
cat(sprintf(
 "- Hostname: %s
  - Start time: %s
  - Populations: %s and %s
  - Output directory: %s
  - Genotyping dataset: %s
  - Gene set annotation file: %s
  - Genome reference file: %s
  - Gene set size limit: %s - %s genes
  - SNP-gene mapping threshold: %skb
  - Permutation cycles: %s\n",
  Sys.info()["nodename"], format(Sys.time(), "%a %b %d %X %Y"),
  popOne, popTwo, popFolder, basename(genotypeFile), basename(annotationFile),
  basename(refgeneFile), minGene, maxGene, snp2geneDist/1000, setPerm
))

######
# (1) RECODE PLINK FAM FILE
######

cat("\n\n(1) ASSIGNING POPULATION STATUS TO GENOTYPE DATA\n\n")
recodeFAM(
  genotypeFile=genotypeFile,
  popOne=popOne,
  popTwo=popTwo,
  populationTable=populationTable,
  outputFile=famName
)
Sys.sleep(3)

######
# (2) DETERMINE POPULATION SUBSTRUCTURE
######

cat("\n\n(2) DETERMINING POPULATION SUBSTRUCTURE\n\n")
popPCA(
  genotypeFile=genotypeFile,
  famFile=famFile,
  popOne=popOne,
  popTwo=popTwo,
  outputFile=pcaFile
)
Sys.sleep(3)

######
# (3) CALCULATE POPULATION-BASED FST
######

cat("\n\n(3) CALCULATING POPULATION-BASED FST\n\n")
calcFST(
  genotypeFile=genotypeFile,
  famFile=famFile,
  popOne=popOne,
  popTwo=popTwo,
  outputFolder=fstFolder,
  outputFile=fstFile
)
Sys.sleep(3)

######
# (4) MAP SNPS TO GENES
######

cat("\n\n(4) MAPPING INPUT SNPS TO NEAREST GENE\n\n")
SNP2gene(
  inFile=snpFile,
  refgeneFile=refgeneFile,
  outputFile=snp2geneFile
)
Sys.sleep(3)

######
# (5) RUN GSEA
######

cat("\n\n(5) RUNNING GENE SET ENRICHMENT ANALYSIS\n\n")
setupGSEArun(
  fstFile=fstFile,
  annotationFile=annotationFile,
  snp2geneFile=snp2geneFile,
  setPerm=setPerm,
  minGene=minGene,
  maxGene=maxGene,
  scoreType=scoreType,
  outputFolder=gseaFolder
)
Sys.sleep(3)

# End sinking
cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%a %b %d %X %Y")))
sink()
