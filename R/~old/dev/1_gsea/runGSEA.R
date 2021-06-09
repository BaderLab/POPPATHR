#' Runs pathway analysis with plink + GSEA using difference in minor allele
#' frequency between two populations (calls calculateMAFdiff.R and
#' setupGSEArun.R)

#-------------------------------------------------------------------------------
## INPUT FILES
pop1 <- "CEU"       # switch based on which pop you want as reference when
pop2 <- "ASW"       # calculating the SNP test statistic (delta MAF)

#dataset <- "HM3"    # choose one or the other
#dataset <- "PNC"
dataset <- "1KG"

geno <- "1kg_phase1_all_autosomes"
fam1 <- "1kg_phase1_all_CEU"
fam2 <- "1kg_phase1_all_ASW"

dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways"
genoF   <- sprintf("%s/data/%s/all/%s", dataDir, dataset, geno)
fam1F <- sprintf("%s/data/%s/CEU/%s.fam", dataDir, dataset, fam1)
fam2F <- sprintf("%s/data/%s/ASW/%s.fam", dataDir, dataset, fam2)

# Path to GWAS .assoc file
assocF <- sprintf(paste0("%s/methods/1_gsea/1KG/CEU_ASW/",
                         "out_170418_pval/gwas/gwas.assoc"), dataDir)

dt	  <- format(Sys.Date(), "%y%m%d")
outDir	<- sprintf("%s/methods/1_gsea/1KG/%s_%s/out_%s_MAFdiff",
                   dataDir, pop1, pop2, dt)

#-------------------------------------------------------------------------------
## SOFTWARE & ANNOTATION
PLINK 		<- "/media/catherine/DATAPART1/Software/plink_linux_1.90/plink"
CALC_GSEA <- "/media/catherine/DATAPART1/Software/GenGen-1.0.1/calculate_gsea.pl"
COMB_GSEA	<- "/media/catherine/DATAPART1/Software/GenGen-1.0.1/combine_gsea.pl"

geneFile <- sprintf("%s/anno/refGene/refGene.hg19.header.txt", dataDir)
pathFile  <- sprintf(paste0("%s/anno/baderlab/Human_GOBP_AllPathways_no_GO_",
                            "iea_April_24_2016_symbol.gmt"), dataDir)

#-------------------------------------------------------------------------------
## WORK BEGINS
if (file.exists(outDir)) {
  cat("Directory exists! Not overwriting\n")
  Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)

#-------------------------------------------------------------------------------
## Determine population substructure
mdsDir <- sprintf("%s/mds", dirname(genoF))
if (!file.exists(mdsDir)) dir.create(mdsDir)
outF <- sprintf("%s/ceu_vs_asw", mdsDir)
ggtitle <- paste("PCA -- complete linkage clustering\n",
                 "European vs. African ancestry populations")

require(ggplot2)
source("popPCA.R")
#popPCA(genoF, fam1F, fam2F, dimensions=3L, outF)

#-------------------------------------------------------------------------------
## Map input SNPs to genes
snpFile <- sprintf("%s/snp2gene.txt", outDir)

require(GWAS2pathway) #package by SP (PatientNetworks)
map_SNP2gene(assocF, geneFile, snpFile)

#-------------------------------------------------------------------------------
## Calculate difference in minor allele frequency (MAF) for SNPs
## between both populations
plinkDir <- sprintf("%s/plink", outDir)
if (!file.exists(plinkDir)) dir.create(plinkDir)

source("calculateMAFdiff.R")
calculateMAFdiff(genoF)

#-------------------------------------------------------------------------------
## Run GSEA
resDir  <- sprintf("%s/pathway_analysis", outDir)
if (!file.exists(resDir)) dir.create(resDir)
sink(sprintf("%s/runGSEA.log", resDir))

# Writes statement to log file indicating whether script is being run
# on personal workstation or Scinet GPC
x <- Sys.info()["nodename"]
if (any(grep("^gpc", x))) cat("on scinet\n") else cat("on workstation\n")

snp2geneF <- sprintf("%s/snp2gene.txt", outDir)
realF   <- sprintf("%s/markerMAF.txt", plinkDir)
statOut <- sprintf("%s/gseaStatFile.txt", resDir)
leOut	  <- sprintf("%s/gseaLEout.txt", resDir)

require(xlsx)
source("setupGSEArun.R")
setupGSEArun(realF, pathFile, snp2geneF,
             statOut=statOut, leOut=leOut,
             CALC_GSEA=CALC_GSEA, COMB_GSEA=COMB_GSEA,
             format4EM=TRUE)

# Specify if you want to format GSEA results for use in Cytoscape
EMdir <- sprintf("%s/enrichment_map", outDir)
if (!file.exists(EMdir)) dir.create(EMdir)

# Close log
sink(NULL)
