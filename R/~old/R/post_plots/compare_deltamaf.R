# [res]
require(ggplot2)
require(reshape2)

## only compare between hm3 and pnc
hm3 <- read.table("out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/freq/markerMAF.txt", h=T, as.is=T)
hm3$dataset <- "HM3"
pnc <- read.table("out_170723_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_MAF_20-200gene/freq/markerMAF.txt", h=T, as.is=T)
pnc$dataset <- "PNC"

both <- rbind(hm3, pnc)

p <- ggplot(both, aes(x=CHI2)) +
        facet_grid(.~dataset, scales="free_x") +
        geom_density(colour="#7fc97f", fill="#7fc97f") +
        ggtitle("Distibution of the SNP-level test statistic") +
        labs(y="Density",
             x=expression(bold(paste(Delta,"MAF")))) +
        theme_Publication()

ggsave("deltamaf_eur-afr_datasets.png", p, width=8, height=4)

### comparison of deltaMAF per SNP from ceu-asw vs. asw-ceu vs. abs
#hm3 files
hm3_ceu_asw <- read.table("out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/freq/markerMAF.txt", h=T, as.is=T)
hm3_asw_ceu <- read.table("out_170927_HM3_pops_hg19_ASW-CEU_MAF_20-200gene/freq/markerMAF.txt", h=T, as.is=T)
hm3_abs <- read.table("out_170724_HM3_pops_hg19_CEU-ASW_absMAF_20-200gene/freq/markerMAF.txt", h=T, as.is=T)

hm3_ceu_asw$def <- "European-African"
hm3_asw_ceu$def <- "African-European"
hm3_abs$def <- "Absolute value"

hm3 <- rbind(hm3_ceu_asw, hm3_asw_ceu, hm3_abs)
hm3_melt <- melt(hm3, id.vars=c("Marker", "def"), measure.vars="CHI2")
hm3_melt$dataset <- "HM3"
#dat2 <- dat_melt$value[dat_melt$value == 0] <- NA
#dat2 <- na.omit(dat_melt)

# pnc files
pnc_ceu_asw <- read.table("out_170723_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_MAF_20-200gene/freq/markerMAF.txt", h=T, as.is=T)
pnc_asw_ceu <- read.table("out_170927_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_MAF_ASW-CEU_20-200gene/freq/markerMAF.txt", h=T, as.is=T)
pnc_abs <- read.table("out_170724_PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW_absMAF_20-200gene/freq/markerMAF.txt", h=T, as.is=T)

pnc_ceu_asw$def <- "European-African"
pnc_asw_ceu$def <- "African-European"
pnc_abs$def <- "Absolute value"

pnc <- rbind(pnc_ceu_asw, pnc_asw_ceu, pnc_abs)
pnc_melt <- melt(pnc, id.vars=c("Marker", "def"), measure.vars="CHI2")
pnc_melt$dataset <- "PNC"
#dat2 <- dat_melt$value[dat_melt$value == 0] <- NA
#dat2 <- na.omit(dat_melt)

both <- rbind(hm3_melt, pnc_melt)

# plot densities of deltaMAF per defintion
p <- ggplot(both, aes(x=value)) +
        facet_grid(dataset~def, scales="free_x") +
        geom_density(colour="#7fc97f", fill="#7fc97f") +
        ggtitle("Distibution of SNP-level test statistic across definitions") +
        labs(y="Density",
             x=expression(bold(paste(Delta,"MAF")))) +
        theme_Publication()

ggsave("compare_deltamaf_defs_datasets.png", p, width=11, height=6)

# MAF per population (MAF_A = cases, MAF_U = controls)
freq <- read.table("out_170724_HM3_pops_hg19_CEU-ASW_MAF_20-200gene/freq/HM3_pops_hg19.frq.cc", h=T, as.is=T)
freq <- na.omit(freq)

maf_ceu <- freq[,c("SNP","MAF_U")]
maf_ceu$population <- "European"
colnames(maf_ceu)[2] <- "MAF"
maf_asw <- freq[,c("SNP","MAF_A")]
maf_asw$population <- "African"
colnames(maf_asw)[2] <- "MAF"

maf_pops <- rbind(maf_ceu, maf_asw)
maf_melt <- melt(maf_pops, id.vars=c("SNP", "population"), measure.vars="MAF")

p <- ggplot(maf_melt, aes(x=value)) +
        facet_grid(.~population) +
        geom_histogram(binwidth=0.03, position="identity", fill="#7fc97f") +
        labs(x="MAF", y="Density") +
        ggtitle("SNP minor allele frequency per population") +
        theme_Publication()

ggsave("ceu_asw_maf.png", p, width=9, height=4)
