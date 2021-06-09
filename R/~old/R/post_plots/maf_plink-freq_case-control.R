#[poppaths]/data
require(ggplot2)

#HM3
#maf_all <- read.table("HM3_2010-05_phase3/HM3_pops_hg19.frq", h=T)
#maf_cc <- read.table("HM3_2010-05_phase3/HM3_pops_hg19_CEU_ASW.frq.cc", h=T)
maf_all <- read.table("HM3/HM3_pops_hg19.frq", h=T)
maf_cc <- read.table("HM3/HM3_pops_hg19_CEU_ASW.frq.cc", h=T)

all <- maf_all[,c("SNP", "MAF")]
all$pop <- "Genome-wide"
ceu <- maf_cc[,c("SNP", "MAF_U")]
colnames(ceu)[2] <- "MAF"
ceu$pop <- "European"
asw <- maf_cc[,c("SNP", "MAF_A")]
colnames(asw)[2] <- "MAF"
asw$pop <- "African"

dat_hm3 <- rbind(all, ceu, asw)
dat_hm3$set <- "HM3"

#PNC
maf_all <- read.table("PNC/all/PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW.frq", h=T)
maf_cc <- read.table("PNC/all/PNC_imputed_merged.CLEAN_FINAL_5sd_noAxiom_CEU_ASW.frq.cc", h=T)

all <- maf_all[,c("SNP", "MAF")]
all$pop <- "Genome-wide"
ceu <- maf_cc[,c("SNP", "MAF_U")]
colnames(ceu)[2] <- "MAF"
ceu$pop <- "European"
asw <- maf_cc[,c("SNP", "MAF_A")]
colnames(asw)[2] <- "MAF"
asw$pop <- "African"

dat_pnc <- rbind(all, ceu, asw)
dat_pnc$set <- "PNC"

# both
both <- rbind(dat_hm3, dat_pnc)
both$MAF <- as.numeric(both$MAF)

# frequency polygraph
title <- paste("Distribution of population-stratified MAF as compared to\n",
               "the genome-wide distribution")
p <- ggplot(both, aes(MAF, colour=pop)) +
      facet_grid(.~set, scales="free_y") +
      geom_freqpoly() +
      ggtitle(title) +
      labs(x="Minor allele frequency", y="Number of SNPs (1000s)") +
      scale_y_continuous(labels=function(x)x/1000) + #divide scale by 1000
      scale_color_manual(values=c("#fdb462", "#984ea3", "#662506")) +
      theme_Publication()

ggsave("maf_dist_genome_pop-strat.png", p, width=8, height=4.5)
