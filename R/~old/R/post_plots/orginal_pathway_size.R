require(reshape2)
require(plyr)
require(dplyr)

source("../../bin/R/post_plots/readPathways.R")
dat <- readPathways("Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt",
    MIN_SIZE=0, MAX_SIZE=1000000)

dat_melt <- melt(dat)
colnames(dat_melt)[2] <- "pathway"

# hc pathways
hc <- read.delim("../../res/out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/ld/hc_snps/pathways_hc.txt", h=T)
hc_path <- as.data.frame(hc[,1])

colnames(hc_path) <- "pathway"
hc_path$pathway <- gsub("\\%.*", "", hc_path$pathway)

path_merge <- join(hc_path, dat_melt)

lengths <- c()
for (i in 1:length(unique(path_merge$pathway))) {
  path <- filter(path_merge, pathway==hc_path[i,])
  lengths[[i]] <- nrow(path)
}

lengths <- as.data.frame(lengths)

final <- data.frame(pathway=hc_path,
                    N_real_pathway=lengths)

write.table(final, "pathways_hc_oringial_path_size.txt", col.names=T,
    row.names=F, quote=F, sep="\t")
