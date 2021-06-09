# only to run LD stats (run in workdir)
dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways"

pop1 <- "CEU"
pop2 <- "YRI"

suppressMessages({
  require(plyr)            #LDstats
  require(dplyr)           #calcMAFdiff, getSNPlists, LDstats, pathPosSel
  require(ggplot2)         #popPCA, LDstats, pathPosSel
  require(GenomicRanges)   #SNP2gene
  require(data.table)      #SNP2gene, calcMAFdiff
  require(tools)           #calcMAFdiff, getPathStats
  require(stringr)         #getPathStats
  require(snpStats)        #LDstats
  require(reshape2)        #LDstats, pathPosSel
  require(gdata)           #LDstats
  require(RColorBrewer)    #LDstats
  require(gridExtra)       #LDstats
  require(cowplot)         #LDstats
  require(splitstackshape) #pathPosSel
  require(utils)           #interactive
  require(SNPRelate)       #LDstats (composite LD)
})

source("../../bin/R/themePublication.R")

ldDir   <- sprintf("%s/ld", getwd())
#hcSnps  <- sprintf("%s/hc_snps", ldDir)
hcSnps  <- sprintf("%s/hc_em_groups", ldDir)
lcSnps <- sprintf("%s/lc_snps", ldDir)

hcInDir=hcSnps
lcInDir=lcSnps
popNames <- c("CEU", "YRI")

statistic <- "R.squared" # or D.prime

# within pathway
statDir <- sprintf("%s/res_%s_%s_updated_script_%s", ldDir,
            basename(hcSnps), basename(lcSnps), statistic)
if (!file.exists(statDir)) dir.create(statDir)

source("../../bin/R/LDstats.R")
LDstats(hcInDir=hcSnps, lcInDir=lcSnps, statistic="R.squared",
        popNames=popNames, outDir=statDir)

# between pathway
statDir <- sprintf("%s/res_%s_%s_BPM_w_singletons", ldDir, basename(hcSnps), basename(lcSnps))
if (!file.exists(statDir)) dir.create(statDir)

snp2geneF <- "gsea/snp2gene.txt"

source("../../bin/R/LDstatsBPM.R")
LDstats(hcInDir=hcSnps, lcInDir=lcSnps, popNames=popNames,
        snp2geneF=snp2geneF, outDir=statDir)
