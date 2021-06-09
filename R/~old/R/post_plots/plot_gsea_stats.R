# [workdir]/gsea/pathway_analysis/
require(reshape2)
require(ggplot2)

allres <- read.delim("results.txt", h=T, as.is=T)
allres$Pathway <- seq.int(nrow(allres))

hc_stats <- read.delim("../../ld/hc_snps/pathways_hc.txt", h=F, as.is=T)
lc_stats <- read.delim("../../ld/lc_snps/pathways_lc.txt", h=F, as.is=T)

allres$pathway_type <- NA
top <- match(hc_stats$V1, allres$Geneset)
bottom <- match(lc_stats$V1, allres$Geneset)

for (i in top) allres[i,]$pathway_type <- "Enriched"
for (i in bottom) allres[i,]$pathway_type <- "Nonenriched"
allres$pathway_type <- as.factor(allres$pathway_type)

blah <- melt(allres, id.vars=c("Pathway", "pathway_type"),
             measure.vars=c("NominalP", "FDR", "ES", "NES"))

p <- ggplot(blah, aes(x=Pathway, y=value, colour=pathway_type)) +
        facet_wrap(~variable, scales = "free_y") +
        geom_point() +
        scale_y_continuous("Statistic value") +
        ggtitle("Comparison of GSEA statistics per enriched and nonenriched pathway") +
        scale_color_discrete(breaks=levels(blah$pathway_type),
                             na.value="lightgrey")
p + theme_Publication()
ggsave("gsea_stat_replicate_high_lowconf.png", width=9, height=8.5)
