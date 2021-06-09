#[workdir]/gsea/pathway_analysis
# Create GSEA enrichment plots using 'fgsea' pkg
# https://bioconductor.org/packages/3.7/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
require(fgsea)
require(ggplot2)
require(gridExtra)
require(dplyr)

# gseaDir
dat <- read.table("../../ld/hc_snps/gseaStat_per_hc_pathway.txt", h=T, as.is=T)

# create vector for gene-level dmaf statistics
vals <- dat[,c(1,3)]
vals <- unique(vals[order(vals$fst, decreasing=T),])
ranks <- setNames(vals$fst, vals$gene)

path_dat <- c()
plot_list = list()
for (i in 1:length(unique(dat$pathway))) {
  path_dat[[i]] <- filter(dat, pathway==i)
  p = plotEnrichment(path_dat[[i]]$gene, ranks) +
        labs(title=sprintf("Pathway_%i", i),
             y="Enrichment score", x="Rank") +
        theme_Publication()
  plot_list[[i]] = p
}

# enrichment plots for all pathways
pdf(file="es_plots_hc.pdf", onefile=T)
for (i in 1:length(unique(dat$pathway))) {
    print(plot_list[[i]])
}
dev.off()

# enrichment plot for hemo and toll pathways
#hemo <- filter(dat, pathway==15)
#toll <- filter(dat, pathway==17)

#p1 <- plotEnrichment(hemo$gene, ranks) +
#        labs(title="Regulation of hemopoiesis",
#             y="Enrichment score", x="Rank") +
#        theme_Publication()

#p2 <- plotEnrichment(toll$gene, ranks) +
#        labs(title="Toll-like receptor signaling pathway",
#             y="Enrichment score", x="Rank") +
#        theme_Publication()

#both <- arrangeGrob(p1, p2, ncol=2, nrow=1)
#ggsave("es_hemo_toll.png", both, width=10, height=4)

# get pathway stats
#hc <- read.delim("../../ld/hc_snps/pathways_hc.txt", h=F)
#hc$Geneset <- gsub(" ", "_", hc$Geneset)
#hc$Geneset <- gsub("\\%.*", "", hc$Geneset)

# vector list of each pathway elements (genes)
#hc_names <- as.data.frame(hc$Geneset)
#hc_names$pathway <- seq.int(nrow(hc_names))

#hc_df <- merge(dat, hc_names, by="pathway")
#hc_df <- hc_df[,c(2,6)]
#hc_df <- unstack(hc_df)

# gsea stats per pathway
#hc_stats <- hc[,c(1,4,5,6)]
#colnames(hc_stats) <- c("pathway", "NES", "pval", "padj")

# table of enrichment plots
#pdf("gseaEStable_hc.pdf", width=13.5, height=10, family='Helvetica')
#plotGseaTable(hc_df, ranks, fgseaRes=hc_stats, gseaParam=0.5)
#dev.off()

###########################
dat <- read.table("../../ld/lc_snps/gseaStat_per_lc_pathway.txt", h=T, as.is=T)

# create vector for gene-level dmaf statistics
vals <- dat[,c(1,3)]
vals <- unique(vals[order(vals$fst, decreasing=T),])
ranks <- setNames(vals$fst, vals$gene)

plot_list = list()
for (i in 1:length(unique(dat$pathway))) {
  path_dat[[i]] <- filter(dat, pathway==i)
  p = plotEnrichment(path_dat[[i]]$gene, ranks) +
        labs(title=sprintf("Pathway_%i", i),
             y="Enrichment score", x="Rank") +
        theme_Publication()
  plot_list[[i]] = p
}

# enrichment plots for all pathways
pdf(file="es_plots_lc.pdf", onefile=T)
for (i in 1:length(unique(dat$pathway))) {
    print(plot_list[[i]])
}
dev.off()

# enrichment plot for enriched pathway 1
#path1 <- filter(dat, pathway==1)
#path1_gene <- path1$gene

# get pathway stats
#lc <- read.delim("../../ld/lc_snps/pathways_lc.txt", h=T)
#lc$Geneset <- gsub(" ", "_", lc$Geneset)
#lc$Geneset <- gsub("\\%.*", "", lc$Geneset)

# vector list of each pathway elements (genes)
#lc_names <- as.data.frame(lc$Geneset)
#lc_names$pathway <- seq.int(nrow(lc_names))

#lc_df <- merge(dat, lc_names, by="pathway")
#lc_df <- lc_df[,c(2,6)]
#lc_df <- unstack(lc_df)

# gsea stats per pathway
#lc_stats <- lc[,c(1,4,5,6)]
#colnames(lc_stats) <- c("pathway", "NES", "pval", "padj")

# table of enrichment plots
#pdf("gseaEStable_lc.pdf", width=13.5, height=10, family='Helvetica')
#plotGseaTable(lc_df, ranks, fgseaRes=lc_stats, gseaParam=0.5)
#dev.off()
