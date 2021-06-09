# [workdir]/ld
require(ggplot2)
require(cowplot)
require(ggpubr)
require(ggsignif)

hc <- read.table("hc_snps/gseaStat_per_hc_pathway.txt", h=T)
lc <- read.table("lc_snps/gseaStat_per_lc_pathway.txt", h=T)
lc$set <- "Unenriched" #instead of 'nonenriched'
#all <- read.table("../gsea/pathway_analysis/gseaStat_all_pathways.txt", h=T, as.is=T)
#all$snp <- as.factor(all$snp)
#all$set <- "All pathways"

dat <- rbind(hc, lc)
dat$pathway <- as.factor(dat$pathway)

p <- ggplot(dat, aes(x=set, y=fst, fill=set)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1, fill="white") +
        labs(y=expression("F"[ST]),
             title="Pathway-level FST\ndistributions") +
        stat_compare_means(method="wilcox", #method "t.test" "wilcox"
                           method.args=list(alternative="less"),
                           label.y=1.05) +
        scale_fill_manual(values=c("#fb9a99", "#386cb0")) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              axis.title.x=element_blank(),
              legend.position="none")

save_plot("fst_dist_enrich_unenrich_wilcox.tiff", p, base_width=3,
    base_height=6.5)

# using ggsignif
pval <- wilcoxon.test(hc$fst, lc$fst)$p.value
anno <- paste(sprintf("p = %s", formatC(pval, digits=2)))
anno2 <- "p < 0.001"

# NOTE changes for publication (2018-12-11)
hc$set <- "Selection-enriched"

dat <- rbind(hc, lc)
dat$pathway <- as.factor(dat$pathway)

dat$set <- factor(dat$set, levels=c("Selection-enriched", "Unenriched"))

p <- ggplot(dat, aes(x=set, y=fst, fill=set)) +
        geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1, fill="white") +
        labs(x=NULL,
             y=expression(paste("Strength of positve selection (F"[ST] ,")"))) +
        geom_signif(annotation=anno2,
                    y_position=1.09, xmin=1, xmax=2,
                    tip_length=0.04,
                    textsize=8) +
        scale_fill_manual(values=c("#fb9a99", "#386cb0")) +
        theme(legend.position="none",
              text=element_text(size=20),
              axis.text.x=element_text(size=20),
              axis.text.y=element_text(size=20))

save_plot("fst_dist_enriched_unenrich_stem_wilcox.svg", p,
    base_width=3.5, base_height=6.5)
