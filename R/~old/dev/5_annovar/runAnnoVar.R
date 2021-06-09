#' Run AnnoVar (calls annotateSnps.R)
#-----------------------------INPUT FILES---------------------------------------
dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways/methods"
snpList <- sprintf(paste0("%s/1_gsea/PNC/CEU_ASW/out_161109_MAFdiff_CEU-ASW/",
                          "post_analyses/snp_list"), dataDir)

dt	  <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/2_annovar/PNC/out_%s_1000g2015aug_all", dataDir, dt)

#------------------------SOFTWARE & ANNOTATION----------------------------------
CONVERT2AV <- "/media/catherine/DATAPART1/Software/annovar/convert2annovar.pl"
TABLE_AV <- "/media/catherine/DATAPART1/Software/annovar/table_annovar.pl"
ANNO_VAR <- "/media/catherine/DATAPART1/Software/annovar/annotate_variation.pl"

annodb <- "/media/catherine/DATAPART1/Software/annovar/humandb_hg19"
dbsnpFile <- sprintf("%s/hg19_avsnp144.txt", annodb)

#--------------------------------WORK BEGINS------------------------------------
if (file.exists(outDir)) {
  cat("Directory exists! Not overwriting\n")
  Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)

# Load additional R libraries
suppressMessages(library(xlsx))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(colorspace))

#-------------------------------RUN ANNOVAR-------------------------------------
source("annotateSnps.R")
annotateSnps(snpList, dbsnpFile,
             snpGroup=c("top","bottom","random","all"),
             TABLE_AV=TABLE_AV, ANNO_VAR=ANNO_VAR,
             annodb=annodb)
