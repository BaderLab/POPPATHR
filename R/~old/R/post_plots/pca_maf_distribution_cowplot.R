# [poppaths]/data/HM3_2010-05_phase3
require(RColorBrewer)
require(ggplot2)
require(cowplot)

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

# minor allele freq plot for ceu and yri
str1 <- sprintf("%s --bed %s.bed --bim %s.bim --fam %s_CEU_YRI.fam", PLINK, dat, dat, dat)
str2 <- sprintf("--allow-no-sex --freq case-control --out %s_CEU_YRI", dat)
cmd <- paste(str1, str2)
system(cmd)

# for ceu and lwk
str1 <- sprintf("%s --bed %s.bed --bim %s.bim --fam %s_CEU_LWK.fam", PLINK, dat, dat, dat)
str2 <- sprintf("--allow-no-sex --freq case-control --out %s_CEU_LWK", dat)
cmd <- paste(str1, str2)
system(cmd)

freq1 <- read.table(sprintf("%s_CEU_YRI.frq.cc", dat), h=T, as.is=T)
freq_ceu <- freq1[,c(2,6)]
colnames(freq_ceu)[2] <- "MAF"
freq_ceu$Population <- "CEU"
freq_yri <- freq1[,c(2,5)]
colnames(freq_yri)[2] <- "MAF"
freq_yri$Population <- "YRI"

freq2 <- read.table(sprintf("%s_CEU_LWK.frq.cc", dat), h=T, as.is=T)
freq_lwk <- freq2[,c(2,5)]
colnames(freq_lwk)[2] <- "MAF"
freq_lwk$Population <- "LWK"

freq_all <- rbind(freq_ceu, freq_yri, freq_lwk)

# maf density plot
p2 <- ggplot(freq_all, aes(x=MAF, colour=Population)) +
        geom_freqpoly() +
        labs(y="Number of SNPs, thousands", x="Minor allele frequency") +
        scale_y_continuous(labels=function(x)x/1000) + #divide scale by 1000
        scale_color_manual(values=c("#7D54A5", "#2D82AF", "#A6CEE3")) +
        theme(legend.position="none")

both <- ggdraw() +
        draw_plot(p2 + theme(legend.justification="bottom"), 0, 0, 1, 1) +
        draw_plot(p1 + theme(legend.justification="top"), 0.43, 0.30, 0.58, 0.67) +
        draw_plot_label(c("A", "B"), c(0, 0.45), c(1, 1), size=15)

# change to png (2018-12-11)
save_plot("pca_maf_distribution_ceu_yri_lwk.png", both,
          base_width=8, base_aspect_ratio=1.2)
# change to pdf (2019-02-21)
save_plot("pca_maf_distribution_ceu_yri_lwk.pdf", both,
          base_width=8, base_aspect_ratio=1.2)
