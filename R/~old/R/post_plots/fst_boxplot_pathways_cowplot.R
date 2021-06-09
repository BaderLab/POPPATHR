#[workdir]
require(stringr)
require(reshape2)
require(ggplot2)
require(cowplot)

# first created gseaStatFile for only hc and lc pathways
hcStatsF <- "ld/hc_snps/gseaStatFile_hc.txt"
no_col <- max(count.fields(hcStatsF))
stats <- readLines(hcStatsF)
stats <- str_split_fixed(stats, "\t", no_col)
snp_stats <- t(stats) #transpose data

hc_fst_list <- list()
for (i in 1:ncol(snp_stats)) {
  path <- snp_stats[-1,i]
  path <- as.data.frame(str_split_fixed(path, ",", 3))
  path[path == ""] <- NA
  path <- na.omit(path)
  hc_fst_list[i] <- as.data.frame(path[,3])
}

hc <- melt(hc_fst_list)
colnames(hc) <- c("fst", "pathway")
hc$pathway <- as.factor(hc$pathway)
#hc$pathway <- paste("pathway", hc$pathway, sep="_")
hc$fst <- as.numeric(as.character(hc$fst))
hc$set <- "Enriched"

########## lc pathways
lcStatsF <- "ld/lc_snps/gseaStatFile_lc.txt"
no_col <- max(count.fields(lcStatsF))
stats <- readLines(lcStatsF)
stats <- str_split_fixed(stats, "\t", no_col)
snp_stats <- t(stats) #transpose data

lc_fst_list <- list()
for (i in 1:ncol(snp_stats)) {
  path <- snp_stats[-1,i]
  path <- as.data.frame(str_split_fixed(path, ",", 3))
  path[path == ""] <- NA
  path <- na.omit(path)
  lc_fst_list[i] <- as.data.frame(path[,3])
}

lc <- melt(lc_fst_list)
colnames(lc) <- c("fst", "pathway")
lc$pathway <- as.factor(lc$pathway)
#lc$pathway <- paste("pathway", lc$pathway, sep="_")
lc$fst <- as.numeric(as.character(lc$fst))
lc$set <- "Unenriched"

both <- rbind(hc, lc)

#FST boxplots
p <- ggplot(both, aes(x=pathway, y=fst)) +
        facet_grid(set ~ ., scales="free_x") +
        geom_boxplot(aes(fill=set)) +
        guides(fill=FALSE) +
        labs(x="Pathway", y="FST") +
        ggtitle(expression(bold(paste("Distribution of SNP-level FST per pathway")))) +
        scale_fill_manual(values=c("#fb9a99", "#386cb0"))

save_plot("fst_per_pathway.tiff", p, base_width=14, base_height=5.5)
