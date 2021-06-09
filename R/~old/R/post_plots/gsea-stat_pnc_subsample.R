#[res]
require(reshape2)
require(ggplot2)
require(plyr)

samp1 <- "out_180101_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_sample_1_MAF_CEU-ASW_20-200gene_sample_1/"
samp2 <- "out_180101_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_sample_2_MAF_CEU-ASW_20-200gene_sample_2/"

# read individual results tables from both datasets
res_samp1 <- read.delim(sprintf("%s/gsea/pathway_analysis/results.txt", samp1),
                        h=T, as.is=T)
res_samp2 <- read.delim(sprintf("%s/gsea/pathway_analysis/results.txt", samp2),
                        h=T, as.is=T)
res_samp1$Dataset <- "sample_1"
res_samp2$Dataset <- "sample_2"

melt_samp1 <- melt(res_samp1, id.vars=c("Geneset", "Dataset"),
                   measure.vars=c("ES", "NES", "NominalP", "FDR"))
melt_samp2 <- melt(res_samp2, id.vars=c("Geneset", "Dataset"),
                   measure.vars=c("ES", "NES", "NominalP", "FDR"))
dat_merge <- merge(melt_samp1, melt_samp2, by=c("Geneset", "variable"))

##########
lm_eqn = function(df) {
    m = lm(value.y ~ value.x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
         list(a = format(coef(m)[1], digits = 2),
              b = format(coef(m)[2], digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}
eq <- ddply(dat_merge,.(variable),lm_eqn)

p <- ggplot(dat_merge, aes(x=value.x, y=value.y)) +
      facet_wrap(~variable, scales="free") +
      geom_smooth(method="lm", se=FALSE, color="red", linetype = "dashed",
                  size=1, formula=y ~ x) +
      geom_point(alpha=0.4) +
      labs(x="PNC sample 1", y="PNC sample 2") +
      ggtitle("GSEA statistics -- PNC subsample") +
      theme_Publication()
ggsave("pnc_subsamp_gsea_concord.png", width=9.5, height=8.5)
