require(qqman)

# SNP CHR BP P
dat <- read.table("HM3_pops_hg19.fst", h=T, as.is=T)
dat <- dat[,c(1:3,5)]
colnames(dat)[3:4] <- c("BP", "P")  # change name of FST col to P
dat <- na.omit(dat)

pdf("fst_manhattan.pdf", width=15)
manhattan(dat, main="Mahattan Plot Fst CEU vs. ASW",
          cex=0.4, cex.axis=0.6, logp=FALSE, ylab="Weir and Cockerham Fst")
dev.off()

x <- as.data.frame(table(dat$CHR))

#[workdir]/ld
#fst <- read.delim("hc_snps/gseaStat_per_hc_pathway.txt", h=T, as.is=T)
fst <- read.delim("lc_snps/gseaStat_per_lc_pathway.txt", h=T, as.is=T)

fst <- fst[,c(1:3)]
colnames(fst) <- c("GENE", "SNP", "P")
snp_pos <- read.delim("../freq/HM3_pops_hg19.fst", h=T, as.is=T)
snp_pos <- snp_pos[,c(1:3)]

dat <- merge(fst, snp_pos, by="SNP")
colnames(dat)[5] <- "BP"
dat <- na.omit(dat)

pdf("fst_manhattan.pdf", width=15)
manhattan(dat, main="Mahattan Plot Fst CEU vs. YRI",
          cex=0.4, cex.axis=0.6, logp=FALSE, ylab="Weir and Cockerham Fst")
dev.off()
