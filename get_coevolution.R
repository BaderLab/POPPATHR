#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# POPPATHR: Population-based pathway analysis of SNP-SNP coevolution
# Special thanks to Shraddha Pai
# Last modified 28 February 2020
#------------------------------------------------------------------------------#

# Loads all packages in a way that allows exporting to child environments
packages <- c("tidyverse", "plyr", "data.table", "reshape2", "gdata", "RColorBrewer",
              "snpStats", "RCy3", "cowplot", "argparse")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Loads functions
Rfun <- list.files(pattern="*.R", path="R", full.names=TRUE)
for (file in Rfun) { source(file) }

######
# PARSES USER ARGUMENTS
######

# Makes argparse object and arguments
parser <- ArgumentParser(description=" POPPATHR: Population-based pathway analysis of SNP-SNP coevolution.")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print verbose output [default TRUE]")
parser$add_argument("-p1", "--population_pair_one", type="character", default="CEU_YRI",
                    help="Names of first population pair tested [default %(default)s]")
parser$add_argument("-p2", "--population_pair_two", type="character", default=NULL,
                    help="Name of second population pair tested; optional [default %(default)s]")
parser$add_argument("-g", "--genotype_file", type="character", default="data/genotypes/HM3_2010_05_phase3",
                    help="Path to PLINK (bed, bim, fam formatted) SNP genotype files [default %(default)s]")
parser$add_argument("-a", "--annotation_file", type="character", default="data/annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                    help="Path to gmt file containing pathway annotations [default %(default)s]")
parser$add_argument("-o", "--output_folder", type="character", default="output",
                    help="Path to output folder [default %(default)s]")
parser$add_argument("--ENRICH_NES", type="integer", default=3,
                    help="NES threshold to define enriched pathway gene sets [default %(default)s]")
parser$add_argument("--UNENRICH_NES", type="integer", default=0.1,
                    help="NES threshold to define unenriched control pathway gene sets [default %(default)s]")
parser$add_argument("--ASSOC_FDR", type="integer", default=0.1,
                    help="FDR value to define significant pathway coevolution [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

## NOTE LINES FOR TESTING
pop_pair_one <- "CEU_YRI"
pop_pair_two <- "CEU_LWK"
genotype_file <- "genotypes/HM3_2010_05_phase3"
annotation_file <- "annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt"
output_folder <- "output"
ENRICH_NES = 3
UNENRICH_NES = 0.1
ASSOC_FDR = 0.05
#######

# Sets important parameters
pop_pair_one <- args$population_pair_one
pop_pair_two <- args$population_pair_two
genotype_file <- args$genotype_file
annotation_file <- args$annotation_file
output_folder <- args$output_folder
ENRICH_NES <- args$ENRICH_NES
UNENRICH_NES <- args$UNENRICH_NES
ASSOC_FDR <- args$ASSOC_FDR

# Define I/O file names based on user input information
pop_one <- unlist(lapply(strsplit(pop_pair_one, "_"), "[", 1))
pop_two <- unlist(lapply(strsplit(pop_pair_one, "_"), "[", 2))
fam_file <- sprintf("%s_%s_%s.fam", genotype_file, pop_one, pop_two)

# Search for both population pairs if specified
if (length(pop_pair_two)) {
  pop_pairs <- paste(c(pop_pair_one, pop_pair_two), collapse="|")
}

# Folders for two population comparisons to validate GSEA results, if specified
results_folders <- list.files(pattern=pop_pairs, path=output_folder, full.names=TRUE)
gsea_folders <- list.files(pattern="gsea", path=results_folders, full.names=TRUE)

# Folders for single population comparison where results will be stored
## NOTE storing results in output folder of population_pair_one analysis (eg. CEU_YRI)
results_folder <- grep(pop_pair_one, results_folders, value=TRUE)
gsea_folder <- grep(pop_pair_one, gsea_folders, value=TRUE)

# Define files from get_enrichment.R
results_file <- list.files(pattern="results.txt", path=gsea_folders, full.names=TRUE)
gseaStat_file <- sprintf("%s/gseaStatFile.txt", gsea_folder)
snp2gene_file <- sprintf("%s/snp2gene.txt", gsea_folder)

# Checks if files exist
if (is.na(pop_one)) {
  stop("Is --population_pair_one (-p1) argument in '[POPULATION_NAME]_[POPULATION_NAME] format?")
}
if (is.na(pop_two)) {
  stop("Is --population_pair_one (-p1) argument in '[POPULATION_NAME]_[POPULATION_NAME] format?")
}
if (unique(is.na(results_file))) {
  stop("Did you run get_enrichment.R? GSEA results were not found!")
}

# Generate subfolders
enrich_folder <- sprintf("%s/enriched", gsea_folder)
if (!file.exists(enrich_folder)) dir.create(enrich_folder)
enrichEM_folder <- sprintf("%s/enrichedEM", gsea_folder)
if (!file.exists(enrichEM_folder)) dir.create(enrichEM_folder)
unenrich_folder <- sprintf("%s/unenriched", gsea_folder)
if (!file.exists(unenrich_folder)) dir.create(unenrich_folder)
assoc_folder <- sprintf("%s/assoc", results_folder)
if (!file.exists(assoc_folder)) dir.create(assoc_folder)
wpm_folder <- sprintf("%s/WPM", assoc_folder)
if (!file.exists(wpm_folder)) dir.create(wpm_folder)
bpm_folder <- sprintf("%s/BPM", assoc_folder)
if (!file.exists(bpm_folder)) dir.create(bpm_folder)

# Pre-define downstream files
EM_group_file <- sprintf("%s/EM_group_list.rda", enrichEM_folder)
EM_file <- sprintf("%s/results_EM.txt", gsea_folder)

#------------------------------------------------------------------------------#
# MESURE SNP-SNP COEVOLUTION WITHIN AND BETWEEN SELECTION-ENRICHED PATHWAYS
#------------------------------------------------------------------------------#

######
# PLOT ENRICHMENT MAP OF SELECTION-ENRICHED PATHWAYS
######

message("\n** PLOTTING ENRICHMENT MAP TO VISUALIZE SELECTION-ENRICHED PATHWAYS **\n")
writeEmapFile(
  results_file=results_file,
  ENRICH_NES=ENRICH_NES,
  out_file=EM_file
)

# NOTE Cytoscape requires absolute paths of input files
EM_file <- normalizePath(EM_file)
annotation_file <- normalizePath(annotation_file)

plotEmap(
  gmt_file=annotation_file,
  EM_file=EM_file,
  EM_group_file=EM_group_file,
  output_folder=enrichEM_folder,
  net_name="selection_enriched",
  image_format="png"
)

######
# WRITING SNP/GENE FILES PER ENRICHED/UNENRICHED PATHWAY
######

message("\n** WRITING SNP/GENE FILES PER SELECTION-ENRICHED/UNENRICHED PATHWAY **\n")
writePathFiles(
  genotype_file=genotype_file,
  results_file=results_file,
  fam_file=fam_file,
  gseaStat_file=gseaStat_file,
  snp2gene_file=snp2gene_file,
  EM_group_file=EM_group_file,
  ENRICH_NES=ENRICH_NES,
  UNENRICH_NES=UNENRICH_NES,
	enrich_folder=enrich_folder,
  enrichEM_folder=enrichEM_folder,
  unenrich_folder=unenrich_folder
)

######
# CALCULATE SNP-SNP COEVOLUTION FOR SELECTION-ENRICHED PATHWAYS
######

message("\n** CALCULATING TRANS-CHROMOSOMAL SNP-SNP CORRELATION **\n")
enrich_folder <- sprintf("%s/pathway_files", enrich_folder)
unenrich_folder <- sprintf("%s/pathway_files", unenrich_folder)

## WPM (within-pathway model)
message("\n** WITHIN-PATHWAY **\n")
SNPassocWPM(
  pop_one=pop_one,
  pop_two=pop_two,
  ASSOC_FDR=ASSOC_FDR,
  enrich_folder=enrich_folder,
  unenrich_folder=unenrich_folder,
  output_folder=wpm_folder
)

## BPM (between-pathway model)
message("\n** BETWEEN-PATHWAY **\n")
SNPassocBPM(
  pop_one=pop_one,
  pop_two=pop_two,
  snp2gene_file=snp2gene_file,
  ASSOC_FDR=ASSOC_FDR,
  enrich_folder=enrichEM_folder,
  unenrich_folder=unenrich_folder,
  output_folder=bpm_folder
)
