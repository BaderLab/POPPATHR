# [poppaths]/data/HM3_2010-05_phase3
require(RColorBrewer)
require(ggplot2)
require(cowplot)
require(plyr)

# pca
PLINK <- "../../../Software/plink_mac/plink"
dat <- "HM3_pops_hg19"

str1 <- sprintf("%s --bfile %s --allow-no-sex", PLINK, dat)
str2 <- sprintf("--cluster --mds-plot %i --out %s", 3L, dat)
cmd <- paste(str1, str2)
system(cmd)

mds <- read.table(sprintf("%s.mds", dat), h=T, as.is=T)
pops <- read.table("relationships_w_pops_041510.txt", h=T, as.is=T)
pops <- pops[,c(2,7)]

mds_pops <- merge(mds, pops, by="IID")
mds_pops$population <- as.factor(mds_pops$population)
colnames(mds_pops)[7] <- "Population"

#change factor levels
mds_pops$Population <- factor(mds_pops$Population,
                              levels=c("YRI", "LWK", "MKK", "ASW",
                                       "CHB", "CHD", "JPT",
                                       "TSI", "CEU", "GIH", "MEX"))

cols <- colorRampPalette(brewer.pal(12, "Paired"))
npal <- cols(length(unique(mds_pops$Population)))

# pca plot
p1 <- ggplot(mds_pops, aes(x=C1, y=C2, colour=Population)) +
        geom_point(shape=4) +
        labs(x="PC1", y="PC2") +
        scale_colour_manual(values=npal)

# filterDist plot for ceu and yri
cyDir <- "../../res/out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/freq"
cyFst <- read.delim(sprintf("%s/HM3_pops_hg19.fst", cyDir), h=T, as.is=T)

# for ceu and lwk
clDir <- "../../res/out_180412_HM3_pops_hg19_CEU-LWK_FST_10-300gene_500kb-dist_10000perm_updated_data/freq"
clFst <- read.delim(sprintf("%s/HM3_pops_hg19.fst", clDir), h=T, as.is=T)

cyFst <- cyFst[,c(2,5)]
cyFst$Population <- "CEU_YRI"

clFst <- clFst[,c(2,5)]
clFst$Population <- "CEU_LWK"

fst_all <- rbind(cyFst, clFst)
fst_all <- na.omit(fst_all)

# get mean fst value per population comparison
mean_dat <- ddply(fst_all, "Population", summarise, fst.mean=mean(FST))

# maf density plot
p2 <- ggplot(fst_all, aes(x=FST, colour=Population)) +
        geom_freqpoly() +
        geom_vline(data=mean_dat, aes(xintercept=fst.mean),
                   linetype="dashed", size=0.6,
                   color=c("#2D82AF", "#A6CEE3")) +
        labs(y="Number of SNPs (1000)", x="FST") +
        scale_y_continuous(labels=function(x)x/1000) + #divide scale by 1000
        scale_color_manual(values=c("#2D82AF", "#A6CEE3")) +
        theme(legend.position="none")

both <- ggdraw() +
        draw_plot(p2 + theme(legend.justification="bottom"), 0, 0, 1, 1) +
        draw_plot(p1 + theme(legend.justification="top"), 0.45, 0.35, 0.55, 0.55) +
        draw_plot_label(c("A", "B"), c(0, 0.45), c(1, 0.92), size=15)

save_plot("pca_fst_distribution_ceu_yri_lwk.tiff", both,
          base_width=8, base_aspect_ratio=1.2)
