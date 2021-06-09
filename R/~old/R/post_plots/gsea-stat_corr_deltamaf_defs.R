#[res]
require(reshape2)
require(plyr)
require(dplyr)
require(ggplot2)

# European-African
ea <- read.delim(paste0("out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/gsea",
                        "/pathway_analysis/results.txt"), h=T, as.is=T)

ea_melt <- melt(ea, id.vars="Geneset",
                measure.vars=c("ES", "NES", "NominalP", "FDR"))
ea_melt$definition <- "European-African"

# African-European
#ae <- read.delim(paste0("out_170927_HM3_pops_hg19_ASW-CEU_MAF_20-200gene/gsea",
#                        "/pathway_analysis/results.txt"), h=T, as.is=T)

#ae_melt <- melt(ae, id.vars="Geneset",
#                measure.vars=c("ES", "NES", "NominalP", "FDR"))
#ae_melt$definition <- "African-European"

# absolute
ab <- read.delim(paste0("out_170724_HM3_pops_hg19_CEU-ASW_absMAF_20-200gene/gsea",
                        "/pathway_analysis/results.txt"), h=T, as.is=T)

ab_melt <- melt(ab, id.vars="Geneset",
                measure.vars=c("ES", "NES", "NominalP", "FDR"))
ab_melt$definition <- "Absolute"

both <- merge(ea_melt, ab_melt, by=c("Geneset", "variable"))

# regression line
lm_eqn = function(df) {
    m = lm(value.y ~ value.x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
         list(a = format(coef(m)[1], digits = 2),
              b = format(coef(m)[2], digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

eq <- ddply(both,.(variable), lm_eqn)
p <- ggplot(both, aes(x=value.x, y=value.y)) +
        facet_wrap(~variable, scales="free") +
        geom_point(alpha=0.3) +
        geom_smooth(method="lm", se=FALSE, color="red", linetype = "dashed", size=1,
                    formula=y~x) +
        labs(x="European-African", y="Absolute") +
        ggtitle(paste("Correlation of GSEA pathway enrichment scores\n",
                      "between SNP statistic definitions")) +
        theme_Publication()

ggsave("gsea_corr_all-stats_dmaf_abs-dmaf.png", p, width=9.5, height=8.5)

## just ES
both_es <- filter(both, variable=="NES")

eq <- ddply(both_es,.(variable), lm_eqn)
p <- ggplot(both_es, aes(x=value.x, y=value.y)) +
        geom_point(alpha=0.3) +
        geom_smooth(method="lm", se=FALSE, color="red", linetype = "dashed", size=1,
                    formula=y~x) +
        geom_text(data=eq, aes(x=-5, y=5, label=V1), parse=TRUE, inherit.aes=FALSE) +
        labs(x="European-African", y="Absolute") +
        ggtitle("Correlation of NES between SNP statistic definitions") +
        theme_Publication()

ggsave("gsea_corr_nes_dmaf_abs-dmaf.png", p, width=6.5, height=6.5)
