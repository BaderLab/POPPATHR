#[workdir]

# concat. all pathway snps to single file
# cat *.bim > [set]_all_pathway_snps.txt

hc <- read.table("ld/hc_snps_noMAF/hc_all_pathway_snps.txt", h=F, as.is=T)
lc <- read.table("ld/lc_snps_noMAF/lc_all_pathway_snps.txt", h=F, as.is=T)
#all <- read.table("ld/all_snps/all_pathway_snps.txt", h=F, as.is=T)

########## rand
#org_dat <- read.table("../../data/HM3/HM3_pops_hg19.bim", h=F, as.is=T)
#org_dat_no_hc <- org_dat[!(org_dat$V2 %in% hc$V2), ]
#set.seed(42)
#samp_hc_size <- org_dat_no_hc[sample(nrow(org_dat_no_hc), nrow(hc), replace=F),]
############

#akey <- read.table("all_plots/akey_selection/akey_selection_regions_722.txt",
#                    row.names=NULL, h=T, sep="\t")
#colnames(akey) <- c("chr", "start", "end", "number", "scans")

#lift hg19 coord to hg19 via UCSC (09/28/2017)
akey <- read.table("all_plots/akey_selection/hglft_genome_5f60_c8d5c0.bed",
                    row.names=NULL, h=F, sep="\t")
colnames(akey) <- c("chr", "start", "end")

akey$start <- as.numeric(gsub(",", "", akey$start))
akey$end <- as.numeric(gsub(",", "", akey$end))

is.between <- function(x, a, b) {
  (x - a) * (b - x) >= 0
}

dat <- lc
true <- is.between(dat$V4, akey$start, akey$end)
nums <- which(true == TRUE)
sel_snps_lc <- dat[nums,]
nrow(sel_snps_lc)

######### convert hg18 to hg19
akey$chr <- paste("chr", akey$chr, sep="")
akey <- akey[,c(1:3)]
write.table(akey, "akey_selection_regions_722.bed", col=F, row=F, quote=F,
            sep="\t", options(scipen=999)) #disable sci notation
