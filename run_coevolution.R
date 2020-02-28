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
mapply(source, Rfun)

######
# PARSES USER ARGUMENTS
######

# Makes argparse object and arguments
parser <- ArgumentParser(description=" POPPATHR: Population-based pathway analysis of SNP-SNP coevolution.")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print verbose output [default TRUE]")
parser$add_argument("-p1", "--pop_one", type="character", default="CEU",
                    help="Population name of first population cohort to test [default %(default)s]")
parser$add_argument("-p2", "--population_two", type="character", default="YRI",
                    help="Population name of second population cohort to test [default %(default)s]")
parser$add_argument("-g", "--genotype_file", type="character", default="genotypes/HM3_2010_05_phase3",
                    help="Path to PLINK (bed, bim, fam formatted) SNP genotype files [default %(default)s]")
parser$add_argument("-p", "--population_table", type="character", default="genotypes/relationships_w_pops_041510.txt",
                    help="Path to table defining population genotypes [default %(default)s]")
parser$add_argument("-a", "--annotation_file", type="character", default="annotations/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
                    help="Path to gmt file containing pathway annotations [default %(default)s]")
parser$add_argument("-r", "--refgene_file", type="character", default="annotations/refGene.hg19.header.txt",
                    help="Path to refGene genome annotation file [default %(default)s]")
parser$add_argument("--SET_PERM", type="integer", default=1000,
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

pop_one <- args$pop_one
pop_two <- args$population_two
genotype_file <- args$genotype_file
population_table <- args$population_table
annotation_file <- args$annotation_file
refgene_file <- args$refgene_file
SET_PERM <- args$SET_PERM
MIN_GENE <- args$MIN_GENE
MAX_GENE <- args$MAX_GENE
SNP2GENE_DIST <- args$SNP2GENE_DIST
ENRICH_NES <- args$ENRICH_NES
UNENRICH_NES <- args$UNENRICH_NES

# Checks if files exist
if (!file.exists(genotype_file)) {
  stop("ERROR: SNP genotype files do not exist!")
}

#------------------------------------------------------------------------------#

## STEP 2: Identify genetic coevolution within and bewteen pathways ##
## Use directory of first population analysis as default (eg. CEU_YRI)
pop_one=pops_1[1]
pop_two=pops_2[1]

## Get GSEA results files from each analysis
### Directories
resDirs  <- list.files(pattern="out", path=getwd(), full.names=TRUE)
gseaDirs <- list.files(pattern="gsea", path=resDirs, full.names=TRUE)
resDir <- resDirs[grep(paste(pop_one, pop_two, sep="_"), resDirs)]
gseaDir <- gseaDirs[grep(paste(pop_one, pop_two, sep="_"), gseaDirs)]
### Files
famF <- sprintf("%s_%s_%s.fam", genotype_file, pop_one, pop_two)
gseaStatF <- sprintf("%s/gseaStatFile.txt", gseaDir, resDir)
snp2geneF <- sprintf("%s/snp2gene.txt", gseaDir, resDir)
resF <- list.files(pattern="results.txt", path=gseaDirs, full.names=TRUE)

## Output directories for inter-chromosomal SNP association analyses
assocDir <- sprintf("%s/assoc", resDir)
enrichDir   <- sprintf("%s/enriched", assocDir)
enrichEmDir <- sprintf("%s/enriched_eMap", assocDir)
unenrichDir <- sprintf("%s/unenriched", assocDir)
if (!file.exists(assocDir)) dir.create(assocDir)
if (!file.exists(enrichDir)) dir.create(enrichDir)
if (!file.exists(enrichEmDir)) dir.create(enrichEmDir)
if (!file.exists(unenrichDir)) dir.create(unenrichDir)

# Plot EnrichmentMap to visualize selection-enriched pathways
message("\n**Plotting EnrichmentMap from GSEA results.\n")
eMapF <- unique(substr(basename(resF), 0, nchar(basename(resF))-4))
eMapF <- sprintf("%s/%s_eMap.txt", gseaDir, eMapF)

writeEmapFile(resF=resF, ENRICH_NES=ENRICH_NES, outF=eMapF)
plotEmap(gmtF=annotation_file, eMapF=eMapF, outDir=enrichEmDir,
         netName="selEnrich", imageFormat="png")

# Write SNP/gene files per enriched and unenriched pathway
message("\n**Generating SNPs lists per enriched and unenriched pathway.\n")
writePathFiles(genotype_file=genotype_file, resF=resF, famF=famF,
               gseaStatF=gseaStatF, snp2geneF=snp2geneF,
               ENRICH_NES=ENRICH_NES, UNENRICH_NES=UNENRICH_NES,
						   enrichDir=enrichDir, enrichEmDir=enrichEmDir,
               unenrichDir=unenrichDir)

# Calculate SNP association statisics for selection-enriched pathways
message("\n**Calculating trans-chromosomal SNP association statistics.\n")

## WPM (within-pathway model)
wpmDir <- sprintf("%s/WPM", assocDir)
if (!file.exists(wpmDir)) dir.create(wpmDir)
SNPassocWPM(enrichDir=enrichDir, unenrichDir=unenrichDir,
            pop_one=pop_one, pop_two=pop_two, outDir=wpmDir)

## BPM (between-pathway model)
bpmDir <- sprintf("%s/BPM", assocDir)
if (!file.exists(bpmDir)) dir.create(bpmDir)
SNPassocBPM(enrichDir=enrichEmDir, unenrichDir=unenrichDir,
            pop_one=pop_one, pop_two=pop_two, snp2geneF=snp2geneF, outDir=bpmDir)

#-------------------------------------------------------------------------------
## STEP 3: Getting gene properties for selection-enriched pathways ##
message("\n**Cross-referencing selection-enriched interactions with BioGRID.\n")
bgridDir <- sprintf("%s/biogrid", outDir)
if (!file.exists(bgridDir)) dir.create(bgridDir)

refgene_file <- sprintf("%s/genes_enrich_unique.txt", enrichDir)
getBiogrid(bgridF=bgridF, refgene_file=refgene_file, outDir=bgridDir)


message("\n**Grabbing selection-enriched gene properties.\n")
propDir <- sprintf("%s/geneProp", outDir)
if (!file.exists(propDir)) dir.create(propDir)

geneDir <- enrichDir
refgene_file  <- sprintf("%s/genes_leadingEdge_unique_hc.txt", geneDir)
clustF <- sprintf("%s/gseaStat_per_hc_pathway_updated.txt", geneDir)
outDir=propDir
