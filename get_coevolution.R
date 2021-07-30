#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# POPPATHR: Population-based pathway analysis of SNP-SNP coevolution
# Special thanks to Shraddha Pai
# Last modified 28 February 2020
#------------------------------------------------------------------------------#

# Load all packages in a way that allows exporting to child environments
packages <- c("tidyverse", "dplyr", "data.table", "reshape2", "gdata",
              "RColorBrewer", "snpStats", "RCy3", "cowplot", "argparse")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Loaddfunctions
Rfun <- list.files(pattern="*.R", path="R", full.names=TRUE)
for (file in Rfun) { source(file) }

######
# PARSE USER ARGUMENTS
######

# Make argparse object and arguments
parser <- ArgumentParser(
    description=paste("POPPATHR: Population-based pathway analysis of SNP-SNP coevolution.",
    "This script determines significant SNP-SNP coevolution within and between pathways.")
)
parser$add_argument("-v", "--verbose", action="store_false",
                    help="Print verbose output")
parser$add_argument("-p1", "--population_pair_one", type="character", default="CEU_YRI",
                    help="Names of first population pair tested [default %(default)s]")
parser$add_argument("-p2", "--population_pair_two", type="character", default=NULL,
                    help="Name of second population pair tested; optional [default %(default)s]")
parser$add_argument("-g", "--genotypeFile", type="character", default="data/genotypes/HM3_2010_05_phase3",
                    help="Path to PLINK (bed, bim, fam formatted) SNP genotype files [default %(default)s]")
parser$add_argument("-a", "--annotationFile", type="character", default="data/annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                    help="Path to gmt file containing pathway annotations [default %(default)s]")
parser$add_argument("-o", "--outputFolder", type="character", default="output",
                    help="Path to output folder [default %(default)s]")
parser$add_argument("--enrichFDR", type="numeric", default=0.05,
                    help="FDR threshold to define enriched pathways [default %(default)s]")
parser$add_argument("--unenrichNES", type="integer", default=0.1,
                    help="NES threshold to define unenriched control pathways [default %(default)s]")
parser$add_argument("--assocFDR", type="integer", default=0.05,
                    help="FDR value to define significant pathway coevolution [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

## NOTE LINES FOR TESTING
popPairOne <- "CEU_YRI"
popPairTwo <- "CEU_LWK"
genotypeFile <- "data/genotypes/HM3_2010_05_phase3"
annotationFile <- "data/annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt"
outputFolder <- "output"
enrichFDR <- 0.05
unenrichNES <- 0.1
assocFDR <- 0.05
#######

# Sets important parameters
popPairOne <- args$population_pair_one
popPairTwo <- args$population_pair_two
genotypeFile <- args$genotypeFile
annotationFile <- args$annotationFile
outputFolder <- args$outputFolder
enrichFDR <- args$enrichFDR
unenrichNES <- args$unenrichNES
assocFDR <- args$assocFDR

# Define I/O file names based on user input information
popOne <- unlist(lapply(strsplit(popPairOne, "_"), "[", 1))
popTwo <- unlist(lapply(strsplit(popPairOne, "_"), "[", 2))
famFile <- sprintf("%s_%s_%s.fam", genotypeFile, popOne, popTwo)

# Search for both population pairs if specified
if (length(popPairTwo)) {
  popPairs <- paste(c(popPairOne, popPairTwo), collapse="|")
}

# Folders for two population comparisons to validate GSEA results, if specified
resultsFolders <- list.files(pattern=popPairs, path=outputFolder, full.names=TRUE)
gseaFolders <- list.files(pattern="gsea", path=resultsFolders, full.names=TRUE)

# Folders for single population comparison where results will be stored
## NOTE storing results in output folder of popPairOne (eg. CEU_YRI)
resultsFolder <- grep(popPairOne, resultsFolders, value=TRUE)
gseaFolder <- grep(popPairOne, gseaFolders, value=TRUE)

# Define files from get_enrichment.R (keep popPairOne/popPairTwo order)
resultsFile <- list.files(pattern="results.txt", path=gseaFolders, full.names=TRUE)
resultsFile <- resultsFile[c(grep(popPairOne, resultsFile), grep(popPairTwo, resultsFile))]
snp2geneFstFile <- sprintf("%s/snp2gene_fstTop.txt", gseaFolder)

# Checks if files exist
if (is.na(popOne)) {
  stop("Is --population_pair_one (-p1) argument in '[POPULATION_NAME]_[POPULATION_NAME] format?")
}
if (is.na(popTwo)) {
  stop("Is --population_pair_one (-p1) argument in '[POPULATION_NAME]_[POPULATION_NAME] format?")
}
if (unique(is.na(resultsFile))) {
  stop("Did you run get_enrichment.R? GSEA results were not found!")
}

# Generate subfolders
enrichFolder <- sprintf("%s/enriched", gseaFolder)
if (!file.exists(enrichFolder)) dir.create(enrichFolder)
enrichEMFolder <- sprintf("%s/enrichedEM", gseaFolder)
if (!file.exists(enrichEMFolder)) dir.create(enrichEMFolder)
unenrichFolder <- sprintf("%s/unenriched", gseaFolder)
if (!file.exists(unenrichFolder)) dir.create(unenrichFolder)
assocFolder <- sprintf("%s/assoc", resultsFolder)
if (!file.exists(assocFolder)) dir.create(assocFolder)
wpmFolder <- sprintf("%s/WPM", assocFolder)
if (!file.exists(wpmFolder)) dir.create(wpmFolder)
bpmFolder <- sprintf("%s/BPM", assocFolder)
if (!file.exists(bpmFolder)) dir.create(bpmFolder)

# Pre-define downstream files
emGroupFile <- sprintf("%s/emGroupList.rda", enrichEMFolder)
emFile <- sprintf("%s/resultsEM.txt", gseaFolder)

#------------------------------------------------------------------------------#
# MESURE SNP-SNP COEVOLUTION WITHIN AND BETWEEN SELECTION-ENRICHED PATHWAYS
#------------------------------------------------------------------------------#

######
# PLOT ENRICHMENT MAP OF SELECTION-ENRICHED PATHWAYS
######

message("\n** PLOTTING ENRICHMENT MAP TO VISUALIZE SELECTION-ENRICHED PATHWAYS **\n")
writeEmapFile(
  resultsFile=resultsFile,
  enrichFDR=enrichFDR,
  outFile=emFile
)

# NOTE Cytoscape requires absolute paths for input files
emFile <- normalizePath(emFile)
annotationFile <- normalizePath(annotationFile)

plotEmap(
  gmtFile=annotationFile,
  emFile=emFile,
  emGroupFile=emGroupFile,
  netName="selection_enriched"
)

######
# WRITING SNP/GENE FILES PER ENRICHED/UNENRICHED PATHWAY
######

message("\n** WRITING SNP/GENE FILES PER SELECTION-ENRICHED/UNENRICHED PATHWAY **\n")
writePathFiles(
  genotypeFile=genotypeFile,
  resultsFile=resultsFile,
  famFile=famFile,
  snp2geneFstFile=snp2geneFstFile,
  emGroupFile=emGroupFile,
  enrichFDR=enrichFDR,
  unenrichNES=unenrichNES,
  enrichFolder=enrichFolder,
  enrichEMFolder=enrichEMFolder,
  unenrichFolder=unenrichFolder
)

######
# CALCULATE SNP-SNP COEVOLUTION FOR SELECTION-ENRICHED PATHWAYS
######

message("\n** CALCULATING TRANS-CHROMOSOMAL SNP-SNP CORRELATION **\n")
enrichEMFolder <- sprintf("%s/pathwayFiles", enrichEMFolder)
unenrichFolder <- sprintf("%s/pathwayFiles", unenrichFolder)

## WPM (within-pathway model)
message("\n** WITHIN-PATHWAY **\n")
SNPassocWPM(
  popOne=popOne,
  popTwo=popTwo,
  snp2geneFile=snp2geneFile,
  assocFDR=assocFDR,
  enrichFolder=enrichEMFolder,
  unenrichFolder=unenrichFolder,
  outputFolder=wpmFolder
)

## BPM (between-pathway model)
message("\n** BETWEEN-PATHWAY **\n")
SNPassocBPM(
  popOne=popOne,
  popTwo=popTwo,
  snp2geneFile=snp2geneFile,
  assocFDR=assocFDR,
  enrichFolder=enrichEMFolder,
  unenrichFolder=unenrichFolder,
  outputFolder=bpmFolder
)
