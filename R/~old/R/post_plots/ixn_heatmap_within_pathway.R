# [workdir]/ld (after running LDstats.R)
require(ggplot2)
require(scales)

dat <- readRDS("res_hc_snps_lc_snps_updated_script_R.squared/hc.diff.pairs.pop2.rds")
dat <- dat[,c(1,4,7,10)]

snpgene <- read.table("../gsea/snp2gene.txt", h=F)
snpgene <- snpgene[,-3]
colnames(snpgene) <- c("snp_1", "gene_1")
dat <- merge(dat, snpgene, by="snp_1")
colnames(snpgene) <- c("snp_2", "gene_2")
dat <- merge(dat, snpgene, by="snp_2")

path_names <- read.delim("hc_snps/pathways_hc.txt", h=F, as.is=T)

plot_list = list()
for (i in 1:length(unique(dat$pathway))) {
    p = ggplot(data=subset(dat, pathway==sprintf("pathway_%i", i)), aes(gene_1, gene_2)) +
          geom_tile(aes(fill=R.squared)) + # background colours are mapped according to the value column
          scale_fill_gradient2(low="white",
                               high=muted("midnightblue")) + # determine the colour
          theme(panel.grid.major.x=element_blank(), #no gridlines
                panel.grid.minor.x=element_blank(),
                panel.grid.major.y=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.background=element_rect(fill="white"), # background=white
                axis.text.x=element_text(angle=90, hjust=1, vjust=1, size=6, face="bold"),
                plot.title=element_text(size=8, face="bold"),
                axis.text.y=element_text(size=6, hjust=1, vjust=1, face="bold")) +
          ggtitle(path_names[i,]) +
          theme(legend.title=element_text(face="bold", size=10)) +
          scale_x_discrete(name="") +
          scale_y_discrete(name="") +
          labs(fill=expression(italic(r^2)))
    plot_list[[i]] = p
}

pdf("wpm_heatmaps_asw.pdf")
for (i in 1:length(unique(dat$pathway))) {
    print(plot_list[[i]])
}
dev.off()
