# res
require(reshape2)
require(plyr)
require(ggplot2)
require(cowplot)

yri <- read.delim("out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/gsea/pathway_analysis/results.txt", h=T, as.is=T)
yri$Populations <- "CEU_YRI"
lwk <- read.delim("out_180412_HM3_pops_hg19_CEU-LWK_FST_10-300gene_500kb-dist_10000perm_updated_data/gsea/pathway_analysis/results.txt", h=T, as.is=T)
lwk$Populations <- "CEU_LWK"

melt_yri <-  melt(yri, id.vars=c("Geneset", "Populations"),
                  measure.vars=c("ES", "NES", "NominalP", "FDR"))
melt_lwk <-  melt(lwk, id.vars=c("Geneset", "Populations"),
                  measure.vars=c("ES", "NES", "NominalP", "FDR"))

dat_merge <- merge(melt_yri, melt_lwk, by=c("Geneset", "variable"))

# regression equation
lm_eqn = function(df) {
    m = lm(value.y ~ value.x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
         list(a = format(coef(m)[1], digits = 2),
              b = format(coef(m)[2], digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}
eq <- ddply(dat_merge,.(variable),lm_eqn)

# plot gsea stat replication
p <- ggplot(dat_merge, aes(x=value.x, y=value.y)) +
        facet_wrap(~variable, scales="free") +
        geom_point(alpha=0.2) +
        geom_smooth(method="lm", se=FALSE, color="red", linetype="dashed", size=1,
                    formula=y~x) +
        labs(x="CEU vs. YRI", y="CEU vs. LWK") +
        ggtitle("GSEA pathway replication statistics - CEU vs. YRI and CEU vs. LWK")

save_plot("gsea_ceu-yri-lwk_replicate.tiff", p, base_width=9.5, base_height=8.5)

# with enriched and unenriched points
hc_stats <- read.delim("out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/ld/hc_snps/pathways_hc.txt", h=F, as.is=T)
lc_stats <- read.delim("out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/ld/lc_snps/pathways_lc.txt", h=F, as.is=T)

yri$pathway_type <- "Other"
top <- match(hc_stats$V1, yri$Geneset)
bottom <- match(lc_stats$V1, yri$Geneset)

for (i in top) yri[i,]$pathway_type <- "Enriched"
for (i in bottom) yri[i,]$pathway_type <- "Unenriched"
yri$pathway_type <- as.factor(yri$pathway_type)
yri$pathway_type <- factor(yri$pathway_type, levels=c("Enriched", "Unenriched",
                                                      "Other"))
lwk$pathway_type <- "Other"
top <- match(hc_stats$V1, lwk$Geneset)
bottom <- match(lc_stats$V1, lwk$Geneset)

for (i in top) lwk[i,]$pathway_type <- "Enriched"
for (i in bottom) lwk[i,]$pathway_type <- "Unenriched"
lwk$pathway_type <- as.factor(lwk$pathway_type)
lwk$pathway_type <- factor(lwk$pathway_type, levels=c("Enriched", "Unenriched",
                                                      "Other"))

########################
melt_yri <-  melt(yri, id.vars=c("Geneset", "Populations", "pathway_type"),
                  measure.vars=c("ES", "NES", "NominalP", "FDR"))
melt_lwk <-  melt(lwk, id.vars=c("Geneset", "Populations", "pathway_type"),
                  measure.vars=c("ES", "NES", "NominalP", "FDR"))

dat_merge <- merge(melt_yri, melt_lwk, by=c("Geneset", "variable", "pathway_type"))

# plot
p <- ggplot(dat_merge, aes(x=value.x, y=value.y, colour=pathway_type)) +
      facet_wrap(~variable, scales="free") +
      geom_smooth(method="lm", se=FALSE, color="black", linetype = "dashed", size=1,
                  formula=y ~ x) +
      geom_point(aes(alpha=pathway_type)) +
      labs(x="CEU vs. YRI", y="CEU vs. LWK") +
      ggtitle("GSEA pathway replication statistics - CEU vs. YRI and CEU vs. LWK") +
      scale_alpha_manual(guide='none', values=list(Enriched=1, Unenriched=1, Other=0.1)) +
      scale_color_manual(values=c("#fb9a99", "#386cb0", "lightgrey"))

save_plot("gsea_ceu-yri-lwk_enrich_nonenrich.tiff", p, base_width=9.5,
  base_height=8.5, base_aspect_ratio=1.2)
