# res
require(reshape2)
require(ggplot2)

# hm3
fst_hm3 <- read.table("out_180320_HM3_pops_hg19_CEU-ASW_FST_20-200gene/freq/markerFST.txt", h=T)
maf_hm3 <- read.table("out_180215_HM3_pops_hg19_CEU-ASW_absMAF_10-300gene_500kb-dist_updated_data/freq/markerMAF.txt", h=T)

colnames(fst_hm3)[2] <- "FST"
colnames(maf_hm3)[2] <- "|dMAF|"
both_hm3 <- merge(fst_hm3, maf_hm3, by="Marker")
both_hm3$set <- "HM3"

#pnc
fst_pnc <- read.table("out_180320_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_FST_CEU-ASW_10-300gene_500kb-dist/freq/markerFST.txt", h=T)
maf_pnc <- read.table("out_180214_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_absMAF_CEU-ASW_10-300gene_500kb-dist/freq/markerMAF.txt", h=T)

colnames(fst_pnc)[2] <- "FST"
colnames(maf_pnc)[2] <- "|dMAF|"
both_pnc <- merge(fst_pnc, maf_pnc, by="Marker")
both_pnc$set <- "PNC"

#both sets
both <- rbind(both_hm3, both_pnc)
dat <- melt(both, id.vars=c("Marker", "set"), measure.vars=c("FST", "|dMAF|"))

title <- expression(bold(paste("Distribution of SNP-level FST vs. absolute ", Delta, "MAF")))
p <- ggplot(dat, aes(value, colour=variable, fill=variable)) +
        facet_grid(. ~ set) +
        geom_density(alpha=0.3) +
        labs(x="Value", y="Density", title=title) +
        scale_fill_manual(values=c("#a6cee3", "#7fc97f")) +
        scale_colour_manual(values=c("#a6cee3", "#7fc97f")) +
        theme_Publication()

ggsave("dist_fst_absDmaf_hm3-pnc.png", p, width=8, height=4.5)
