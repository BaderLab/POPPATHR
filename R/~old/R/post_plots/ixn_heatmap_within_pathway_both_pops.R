# [workdir]/ld (after running LDstats.R)
require(ggplot2)
require(scales)

ceu <- readRDS("res_hc_snps_lc_snps_updated_script_R.squared/hc.diff.pairs.pop1.rds")
yri <- readRDS("res_hc_snps_lc_snps_updated_script_R.squared/hc.diff.pairs.pop2.rds")

dat <- rbind(ceu, yri)
dat <- na.omit(dat)

snpgene <- read.table("../gsea/snp2gene.txt", h=F)
snpgene <- snpgene[,-3]
colnames(snpgene) <- c("snp_1", "gene_1")
dat <- merge(dat, snpgene, by="snp_1")
colnames(snpgene) <- c("snp_2", "gene_2")
dat <- merge(dat, snpgene, by="snp_2")

path_names <- read.delim("hc_snps/pathways_hc.txt", h=F, as.is=T)

# replace all punctuation with underscores
path_names$V1 <- gsub('([[:punct:]])|\\s+', '_', path_names$V1)

plot_list = list()
for (i in 1:length(unique(dat$pathway))) {
    p = ggplot(data=subset(dat, pathway==sprintf("pathway_%i", i)), aes(gene_1, gene_2)) +
          facet_grid(. ~ pop) +
          geom_tile(aes(fill=R.squared)) + # background colours are mapped according to the value column
          scale_fill_gradient2(low="white",
                               high=muted("midnightblue")) + # determine the colour
          theme_linedraw() +
          theme(panel.background=element_rect(fill="white"), # background=white
                panel.grid=element_blank(),
                plot.title=element_text(size=12, face="bold", hjust=0.5),
                axis.text.x=element_text(angle=90, hjust=1, vjust=1, size=6, face="bold"),
                axis.text.y=element_text(size=6, hjust=1, vjust=1, face="bold"),
                legend.title=element_text(face="bold", size=10)) +
          ggtitle(path_names[i,]) +
          scale_x_discrete(name="") +
          scale_y_discrete(name="") +
          labs(fill=expression(italic(r^2)))
    plot_list[[i]] = p
}

pdf("wpm_enriched_heatmaps_ceu_yri.pdf", height=13, width=22)
for (i in 1:length(unique(dat$pathway))) {
    print(plot_list[[i]])
}
dev.off()
