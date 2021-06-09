# Run Fst and R2 calculations (calls calculateFst.R / calculateR2.R)

#-----------------------------INPUT FILES---------------------------------------
dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways"

plinkPop1F <- sprintf(paste0("%s/data/PNC/CEU/PNC_imputed_merged.CLEAN_FINAL",
                             "_5sd_noAxiom_CEU"), dataDir)

plinkPop2F <- sprintf(paste0("%s/data/PNC/ASW/PNC_imputed_merged.CLEAN_FINAL",
                             "_5sd_noAxiom_ASW"), dataDir)

hcPop1Dir <- sprintf("%s/input_snps/high_conf/pnc/ceu", dataDir)
hcPop2Dir <- sprintf("%s/input_snps/high_conf/pnc/asw", dataDir)

PLINK 		<- "/media/catherine/DATAPART1/Software/plink_linux_1.90/plink"

dt	  <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/methods/2_snpStats/out_%s_pnc", dataDir, dt)

#--------------------------WORK BEGINS------------------------------------------
suppressMessages({
  require(snpStats); require(plyr); require(dplyr); require(tools)
  require(reshape2); require(ggplot2); require(scales); require(xlsx)
})

if (file.exists(outDir)) {
  cat("Directory exists! Not overwriting\n")
  Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)

#----------------------------RUN FUNCTIONS--------------------------------------
source("calculateAssoc.R")
calculateAssoc(plinkF, highSnpDir, makePlots=FALSE)
