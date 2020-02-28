#!/usr/bin/env Rscript

#------------------------------------------------------------------------------#
# POPPATHR: Population-based pathway analysis of SNP-SNP coevolution
# Special thanks to Shraddha Pai
# Last modified 16 October 2019
#------------------------------------------------------------------------------#

# Loads all packages in a way that allows exporting to child environments
packages <- c("tidyverse", "data.table", "reshape2", "gdata", "RColorBrewer",
              "gridExtra", "cowplot", "GenomicRanges", "snpStats", "RCy3", "argparse")
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
parser$add_argument("-p1", "--population_one", type="character", default="CEU",
                    help="Population name of first population cohort to test [default %(default)s]")
parser$add_argument("-p2", "--population_two", type="character", default="YRI",
                    help="Population name of second population cohort to test [default %(default)s]")
parser$add_argument("-g", "--genotype_file", type="character", default="genotypes/HM3_2010_05_phase3",
                    help="Path to PLINK (bed, bim, fam formatted) SNP genotype files [default %(default)s]")
parser$add_argument("-t", "--population_table", type="character", default="genotypes/relationships_w_pops_041510.txt",
                    help="Path to table defining population genotypes [default %(default)s]")
parser$add_argument("-a", "--annotation_file", type="character", default="annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                    help="Path to gmt file containing pathway annotations [default %(default)s]")
parser$add_argument("-r", "--refgene_file", type="character", default="annotations/refGene.hg19.header.txt",
                    help="Path to refGene genome annotation file [default %(default)s]")
parser$add_argument("-o", "--output_folder", type="character", default="output",
                    help="Path to output folder [default %(default)s]")
parser$add_argument("--SET_PERM", type="integer", default=10000,
                    help="Number of GSEA permutation cycles to run [default %(default)s]")
parser$add_argument("--MIN_GENE", type="integer", default=10,
                    help="Minimum number of genes permitted in pathway gene set [default %(default)s]")
parser$add_argument("--MAX_GENE", type="integer", default=300,
                    help="Maximum number of genes permitted in pathway gene set [default %(default)s]")
parser$add_argument("--SNP2GENE_DIST", type="integer", default=500000,
                    help="Maximum distance (bp) considered for SNP-to-gene mapping [default %(default)s]")
parser$add_argument("--ENRICH_NES", type="integer", default=3,
                    help="NES threshold to define enriched pathway gene sets [default %(default)s]")
parser$add_argument("--UNENRICH_NES", type="integer", default=0.1,
                    help="NES threshold to define unenriched control pathway gene sets [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

pop_one <- args$population_one
pop_two <- args$population_two
genotype_file <- args$genotype_file
population_table <- args$population_table
annotation_file <- args$annotation_file
refgene_file <- args$refgene_file
output_folder <- args$output_folder
SET_PERM <- args$SET_PERM
MIN_GENE <- args$MIN_GENE
MAX_GENE <- args$MAX_GENE
SNP2GENE_DIST <- args$SNP2GENE_DIST
ENRICH_NES <- args$ENRICH_NES
UNENRICH_NES <- args$UNENRICH_NES

# Checks if files exist
snp_file <- sprintf("%s.bim", genotype_file)

if (!file.exists(snp_file)) {
  stop("ERROR: PLINK-formatted SNP genotype file does not exist!")
}

# Generate parent output folder if it doesn't already exist
if (!dir.exists(output_folder)) { dir.create(output_folder) }

# Generate population-specific subfolders
pop_folder <- sprintf("%s/%s_%s", output_folder, pop_one, pop_two)
if (!dir.exists(pop_folder)) { dir.create(pop_folder) }
pca_folder <- sprintf("%s/pca", pop_folder)
if (!file.exists(pca_folder)) dir.create(pca_folder)
fst_folder <- sprintf("%s/fst", pop_folder)
if (!file.exists(fst_folder)) dir.create(fst_folder)
gsea_folder <- sprintf("%s/gsea", pop_folder)
if (!file.exists(gsea_folder)) dir.create(gsea_folder)

#------------------------------------------------------------------------------#
# IDENTIFY SELECTION-ENRICHED PATHWAY GENE SETS
#------------------------------------------------------------------------------#

######
# WORK BEGINS
######

# Sink messages to file
#sink(sprintf("%s/runPipeline.log", pop_folder))

# List user parameters
cat("Running POPPATHR with the following input data and parameters...\n")
cat(sprintf(
 "\tPopulations: %s and %s
  \tOutput directory: %s
  \tGenotyping dataset: %s
  \tPathway size limit: %s - %s genes
  \tSNP-gene mapping threshold: %skb
  \tPermutation cycles: %s
  \tEnriched gene-set NES threshold: %s
  \tUnenriched control gene-set NES threshold: %s\n",
  pop_one, pop_two, pop_folder, genotype_file, MIN_GENE, MAX_GENE,
  SNP2GENE_DIST/1000, SET_PERM, ENRICH_NES, UNENRICH_NES
))
Sys.sleep(3)

# Define I/O file names based on user input information
fam_name <- sprintf("%s_%s_%s", basename(genotype_file), pop_one, pop_two)
fam_file <- sprintf("%s_%s_%s.fam", genotype_file, pop_one, pop_two)
pca_file <- sprintf("%s/%s", pca_folder, fam_name)
fst_file <- sprintf("%s/markerFST.txt", fst_folder)
snp2gene_file <- sprintf("%s/snp2gene.txt", gsea_folder)

######
# RECODE PLINK FAM FILE
######

message("\n** ASSIGNING POPULATION STATUS TO GENOTYPE DATA **\n")
recodeFAM(
  genotype_file=genotype_file,
  pop_one=pop_one,
  pop_two=pop_two,
  population_table=population_table,
  out_file=fam_name
)
Sys.sleep(3)

######
# DETERMINE POPULATION SUBSTRUCTURE
######

message("\n** DETERMINING POPULATION SUBSTRUCTURE **\n")
popPCA(
  genotype_file=genotype_file,
  fam_file=fam_file,
  pop_one=pop_one,
  pop_two=pop_two,
  out_file=pca_file
)
Sys.sleep(3)

######
# CALCULATE POPULATION-BASED FST
######

message("\n** CALCULATING POPULATION-BASED FST **\n")
calcFST(
  genotype_file=genotype_file,
  fam_file=fam_file,
  output_folder=fst_folder,
  out_file=fst_file
)
Sys.sleep(3)

######
# MAP SNPS TO GENES
######

message("\n** MAPPING INPUT SNPS TO NEAREST GENE **\n")
SNP2gene(
  in_file=snp_file,
  refgene_file=refgene_file,
  out_file=snp2gene_file
)
Sys.sleep(3)

######
# RUN GSEA
######

message("\n** RUNNING GENE-SET ENRICHMENT ANALYSIS **\n")
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
