fam <- read.table("HM3_pops_hg19.fam")
pops <- read.delim("relationships_w_pops_041510.txt", h=T)

ceu <- pops[,2][pops[,7] == "CEU"]
ceu1 <- ceu[1:83]
ceu2 <- ceu[84:165]

asw <- pops[,2][pops[,7] == "ASW"]
asw1 <- ceu[1:44]
asw2 <- ceu[45:87]

fam_ceu <- fam
fam_ceu[fam_ceu$V2 %in% ceu1, 'V6'] <- 1
fam_ceu[fam_ceu$V2 %in% ceu2, 'V6'] <- 2

fam_asw <- fam
fam_asw[fam_asw$V2 %in% asw1, 'V6'] <- 1
fam_asw[fam_asw$V2 %in% asw2, 'V6'] <- 2

write.table(fam_ceu, "HM3_pops_hg19_CEU_CEU.fam", col=F, row=F, quote=F)
write.table(fam_asw, "HM3_pops_hg19_ASW_ASW.fam", col=F, row=F, quote=F)


# fst calculations separately for each
require(data.table)
require(reshape2)
require(ggplot2)

ceu <- fread("HM3_pops_hg19_CEU_CEU.fst", h=T, data.table=F)
asw <- fread("HM3_pops_hg19_ASW_ASW.fst", h=T, data.table=F)
ceu_asw <- fread("HM3_pops_hg19_CEU_ASW.fst", h=T, data.table=F)
asw_ceu <- fread("HM3_pops_hg19_ASW_CEU.fst", h=T, data.table=F)

ceu <- ceu[,c(2,5)]
ceu$pop <- "CEU"
asw <- asw[,c(2,5)]
asw$pop <- "ASW"
ceu_asw <- ceu_asw[,c(2,5)]
ceu_asw$pop <- "CEU_ASW"
asw_ceu <- asw_ceu[,c(2,5)]
asw_ceu$pop <- "ASW_CEU"

all <- rbind(ceu, asw, ceu_asw, asw_ceu)
dat <- melt(all, id.vars=c("SNP", "pop"), measure.vars="FST")

p <- ggplot(dat, aes(value, fill=pop)) +
      geom_density(alpha=0.4) +
      labs(x="FST", title="FST CEU vs. ASW pop comparison") +
      xlim(0, 0.4) +
      theme_Publication()

ggsave("fst_pop_compare.png", p, width=7, height=5)
