# [workdir]/gsea/pathway_analysis
# plot pathway size vs. enrichment score
require(reshape2)
require(ggplot2)

res <- read.delim("results.txt", h=T, as.is=T)
hc <- read.delim("pathways_hc.txt", h=F, col.names="Geneset", as.is=T)
lc <- read.delim("pathways_lc.txt", h=F, col.names="Geneset", as.is=T)

res$pathway_type <- "Other"
top <- match(hc$Geneset, res$Geneset)
bottom <- match(lc$Geneset, res$Geneset)

for (i in top) res[i,]$pathway_type <- "Enriched"
for (i in bottom) res[i,]$pathway_type <- "Nonenriched"
res$pathway_type <- as.factor(res$pathway_type)

res_melt <-  melt(res, id.vars=c("Geneset", "Size", "pathway_type"),
                  measure.vars=c("ES", "NES"))

p <- ggplot(res_melt, aes(x=Size, y=value, colour=pathway_type)) +
      facet_wrap(~variable, scales="free") +
      geom_point(aes(alpha=pathway_type)) +
      labs(x="Pathway size", y="Pathway enrichment statistic") +
      ggtitle("Correlation of pathway size vs. enrichment scores") +
      scale_alpha_manual(guide='none', values=list(Enriched=1, Nonenriched=1, Other=0.1)) +
      scale_color_manual(values=c("#fb9a99", "#386cb0", "lightgrey")) +
      theme_Publication()
ggsave("pathway_size_vs_ES_NES.png", p, width=9.5, height=4.5)


####
# [res]
dat <- read.delim("out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/gsea/pathway_analysis/results.txt", h=T)
sel <- read.delim("out_180130_HM3_pops_hg19_CEU-ASW_MAF_20-200gene_selection_genes/gsea/pathway_analysis/results.txt", h=T)

dat$Geneset <- gsub("\\%.*", "", dat$Geneset)
merged <- merge(dat, sel, by="Geneset")
