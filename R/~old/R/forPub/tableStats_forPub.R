library(plyr)
library(dplyr)
dataDir <- "~/PopulationPathways/res/out_180412_HM3_pops_hg19_CEU-YRI_FST_10-300gene_500kb-dist_10000perm_updated_data/"

## WITHIN-PATHWAY
setwd(sprintf("%s/ld/res_hc_snps_lc_snps_updated_script_R.squared/", dataDir))
hc_ceu <- readRDS("hc.diff.pairs.pop1.rds")
hc_yri <- readRDS("hc.diff.pairs.pop2.rds")
lc_ceu <- readRDS("lc.diff.pairs.pop1.rds")
lc_yri <- readRDS("lc.diff.pairs.pop2.rds")

# unenriched stats (easier)
dat <- lc_ceu
dat <- lc_yri

nrow(dat) # num ixns
mean(dat$R.squared) # mean r2
max(dat$R.squared) # top r2
quantile(dat$R.squared, .99) # 99th pctl

# enriched stats
dat <- hc_ceu
dat <- hc_yri

# group by pathway
hc_list <- list()
for (i in 1:length(unique(dat$pathway))) {
 hc_list[[i]] <- dat[dat$pathway==sprintf("pathway_%i", i),]
}

num <- lapply(hc_list, nrow) # num ixns
meanr <- lapply(hc_list, function(x) mean(x$R.squared)) # mean r2
topr <- lapply(hc_list, function(x) max(x$R.squared)) # top r2
prct <- lapply(hc_list, function(x) quantile(x$R.squared, .99)) # 99th pctl

ind <- c(16,14,27:29,7,15) #ceu
ind <- c(13,52,53,17,55,54,7,28,29,23,18,21,47,15,3,45,30,43,48,40,37,
         22,51,4,27,1,26,9,44,46,6,10,20,5,19,35) #yri

print(as.data.frame(unlist(num[ind])),row.names=F)
print(as.data.frame(unlist(meanr[ind])),row.names=F)
print(as.data.frame(unlist(topr[ind])),row.names=F)
print(as.data.frame(unlist(prct[ind])),row.names=F)

## BETWEEN-PATHWAY
setwd(sprintf("%s/ld/res_hc_em_groups_lc_snps_BPM_w_singletons/", dataDir))
hc_ceu <- readRDS("hc.diff.pairs.pop1.rds")
hc_yri <- readRDS("hc.diff.pairs.pop2.rds")
lc_ceu <- readRDS("lc.diff.pairs.pop1.rds")
lc_yri <- readRDS("lc.diff.pairs.pop2.rds")

# unenriched stats
dat <- lc_ceu
dat <- lc_yri

nrow(dat) # num ixns
mean(dat$R.squared, na.rm=T) # mean r2
max(dat$R.squared, na.rm=T) # top r2
quantile(dat$R.squared, .99, na.rm=T) # 99th pctl

# enriched stats
hc_ceu_sig <- read.delim("hc_pvals_per_interaction_alt-l_CEU.txt", h=T, as.is=T)
colnames(hc_ceu_sig)[1:2] <- c("pathway_pair1", "pathway_pair2")
hc_yri_sig <- read.delim("hc_pvals_per_interaction_alt-l_YRI.txt", h=T, as.is=T)
colnames(hc_yri_sig)[1:2] <- c("pathway_pair1", "pathway_pair2")

hc_ceu2 <- join(hc_ceu, hc_ceu_sig, by=c("pathway_pair1", "pathway_pair2"))
hc_ceu2 <- filter(hc_ceu2, fdr <= 0.2)
hc_ceu2 <- hc_ceu2[order(hc_ceu2$fdr, hc_ceu2$pvals),]
hc_yri2 <- join(hc_yri, hc_yri_sig, by=c("pathway_pair1", "pathway_pair2"))
hc_yri2 <- filter(hc_yri2, fdr <= 0.2)
hc_yri2 <- hc_yri2[order(hc_yri2$fdr, hc_yri2$pvals),]

dat <- hc_ceu2
dat <- hc_yri2

# group by pathway
hc_list <- list()
for (i in unique(dat$ixn_num)) {
 hc_list[[i]] <- dat[dat$ixn_num==i,]
}

num <- lapply(hc_list, nrow) # num ixns
meanr <- lapply(hc_list, function(x) mean(x$R.squared, na.rm=T)) # mean r2
topr <- lapply(hc_list, function(x) max(x$R.squared, na.rm=T)) # top r2
prct <- lapply(hc_list, function(x) quantile(x$R.squared, .99, na.rm=T)) # 99th pctl

print(as.data.frame(unlist(num)),row.names=F)
print(as.data.frame(unlist(meanr)),row.names=F)
print(as.data.frame(unlist(topr)),row.names=F)
print(as.data.frame(unlist(prct)),row.names=F)
