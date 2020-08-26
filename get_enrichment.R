#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# POPPATHR: Population-based pathway analysis of SNP-SNP coevolution
# Special thanks to Shraddha Pai
# Last modified 28 February 2020
#------------------------------------------------------------------------------#

# Loads all packages in a way that allows exporting to child environments
packages <- c("tidyverse", "data.table", "reshape2", "gdata", "RColorBrewer",
              "GenomicRanges", "argparse")
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
parser <- ArgumentParser(description=paste("POPPATHR: Population-based pathway analysis of SNP-SNP coevolution.",
                                           "This script identifies selection-enriched pathways between two population cohorts."))
parser$add_argument("-v", "--verbose", action="store_false",
                    help="Print verbose output")
parser$add_argument("-p", "--population_pair", type="character",
                    help="Names of two population cohorts to test (e.g., CEU_YRI or CEU_LWK)")
parser$add_argument("-g", "--genotype_file", type="character", default="data/genotypes/HM3_2010_05_phase3",
                    help="Path to PLINK (bed, bim, fam formatted) SNP genotype files [default %(default)s]")
parser$add_argument("-t", "--population_table", type="character", default="data/genotypes/relationships_w_pops_041510.txt",
                    help="Path to table defining population genotypes [default %(default)s]")
parser$add_argument("-a", "--annotation_file", type="character", default="data/annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                    help="Path to gmt file containing gene set annotations [default %(default)s]")
parser$add_argument("-r", "--refgene_file", type="character", default="data/annotations/refGene.hg19.header.txt",
                    help="Path to refGene genome annotation file [default %(default)s]")
parser$add_argument("-o", "--output_folder", type="character", default="output",
                    help="Path to output folder [default %(default)s]")
parser$add_argument("--SET_PERM", type="integer", default=10000,
                    help="Number of GSEA permutation cycles to run [default %(default)s]")
parser$add_argument("--MIN_GENE", type="integer", default=10,
                    help="Minimum number of genes permitted in gene set [default %(default)s]")
parser$add_argument("--MAX_GENE", type="integer", default=300,
                    help="Maximum number of genes permitted in gene set [default %(default)s]")
parser$add_argument("--SNP2GENE_DIST", type="integer", default=500000,
                    help="Maximum distance (bp) considered for SNP-to-gene mapping [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

# Sets important parameters
pop_pair <- args$population_pair
genotype_file <- args$genotype_file
population_table <- args$population_table
annotation_file <- args$annotation_file
refgene_file <- args$refgene_file
output_folder <- args$output_folder
SET_PERM <- args$SET_PERM
MIN_GENE <- args$MIN_GENE
MAX_GENE <- args$MAX_GENE
SNP2GENE_DIST <- args$SNP2GENE_DIST

# Define I/O file names based on user input information
pop_one <- unlist(lapply(strsplit(pop_pair, "_"), "[", 1))
pop_two <- unlist(lapply(strsplit(pop_pair, "_"), "[", 2))
snp_file <- sprintf("%s.bim", genotype_file)
fam_name <- sprintf("%s_%s", basename(genotype_file), pop_pair)
fam_file <- sprintf("%s_%s.fam", genotype_file, pop_pair)

# Checks if files exist
if (!file.exists(snp_file)) {
  stop("PLINK-formatted SNP genotype file does not exist!")
}
if (is.na(pop_one)) {
  stop("Is --population_pair (-p) argument in '[POPULATION_NAME]_[POPULATION_NAME]' format?")
}
if (is.na(pop_two)) {
  stop("Is --population_pair (-p) argument in '[POPULATION_NAME]_[POPULATION_NAME]' format?")
}

# Create parent output folder if it doesn't already exist
if (!dir.exists(output_folder)) { dir.create(output_folder) }

# Create analysis-specific subfolders
pop_folder <- sprintf("%s/%s", output_folder, pop_pair)
if (!dir.exists(pop_folder)) { dir.create(pop_folder) }
pca_folder <- sprintf("%s/pca", pop_folder)
if (!file.exists(pca_folder)) dir.create(pca_folder)
fst_folder <- sprintf("%s/fst", pop_folder)
if (!file.exists(fst_folder)) dir.create(fst_folder)
gsea_folder <- sprintf("%s/gsea", pop_folder)
if (!file.exists(gsea_folder)) dir.create(gsea_folder)

# Pre-define downstream files
pca_file <- sprintf("%s/%s", pca_folder, fam_name)
fst_file <- sprintf("%s/markerFST.txt", fst_folder)
snp2gene_file <- sprintf("%s/snp2gene.txt", gsea_folder)

#------------------------------------------------------------------------------#
# IDENTIFY SELECTION-ENRICHED GENE SETS
#------------------------------------------------------------------------------#

######
# WORK BEGINS
######

# Sink cats to file
sink(sprintf("%s/get_enrichment.log", pop_folder))

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
  pop_one, pop_two, pop_folder, basename(genotype_file), basename(annotation_file),
  basename(refgene_file), MIN_GENE, MAX_GENE, SNP2GENE_DIST/1000, SET_PERM
))

######
# (1) RECODE PLINK FAM FILE
######

cat("\n\n(1) ASSIGNING POPULATION STATUS TO GENOTYPE DATA\n\n")
recodeFAM(
  genotype_file=genotype_file,
  pop_one=pop_one,
  pop_two=pop_two,
  population_table=population_table,
  output_file=fam_name
)
Sys.sleep(3)

######
# (2) DETERMINE POPULATION SUBSTRUCTURE
######

cat("\n\n(2) DETERMINING POPULATION SUBSTRUCTURE\n\n")
popPCA(
  genotype_file=genotype_file,
  fam_file=fam_file,
  pop_one=pop_one,
  pop_two=pop_two,
  output_file=pca_file
)
Sys.sleep(3)

######
# (3) CALCULATE POPULATION-BASED FST
######

cat("\n\n(3) CALCULATING POPULATION-BASED FST\n\n")
calcFST(
  genotype_file=genotype_file,
  fam_file=fam_file,
  output_folder=fst_folder,
  output_file=fst_file
)
Sys.sleep(3)

######
# (4) MAP SNPS TO GENES
######

cat("\n\n(4) MAPPING INPUT SNPS TO NEAREST GENE\n\n")
SNP2gene(
  in_file=snp_file,
  refgene_file=refgene_file,
  output_file=snp2gene_file
)
Sys.sleep(3)

######
# (5) RUN GSEA
######

cat("\n\n(5) RUNNING GENE SET ENRICHMENT ANALYSIS\n\n")
setupGSEArun(
  fst_file=fst_file,
  annotation_file=annotation_file,
  snp2gene_file=snp2gene_file,
  SNP2GENE_DIST=SNP2GENE_DIST,
  MIN_GENE=MIN_GENE,
  MAX_GENE=MAX_GENE,
  SET_PERM=SET_PERM,
  output_folder=gsea_folder
)
Sys.sleep(3)

# End sinking
cat(sprintf("\nEnd time: %s\n", format(Sys.time(), "%a %b %d %X %Y")))
sink()
