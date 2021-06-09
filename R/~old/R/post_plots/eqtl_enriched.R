#wget http://www.exsnp.org/data/GSexSNP_allc_allp_ld8.txt
require(dplyr)

# High-confidence eQTLs considering all 16 studies (LD > 0.8):
eqtl <- read.delim("GSexSNP_allc_allp_ld8.txt", h=T, as.is=T)
colnames(eqtl)[3] <- "exGene"

# Disease associated High-confidence (r2>0.8) eQTLs :
# High-confidence eQTLs at 0.05 cM distance to disease risk SNPs
hypertension <- read.delim("WTCCC_DZ_allc_allp_GSeQTL_0.05cM_list_HT.txt", h=T, as.is=T)
hypertension$Disease <- "Hypertension"
bipolar <- read.delim("WTCCC_DZ_allc_allp_GSeQTL_0.05cM_list_BD.txt", h=T, as.is=T)
bipolar$Disease <- "Bipolar disorder"
rarthritis <- read.delim("WTCCC_DZ_allc_allp_GSeQTL_0.05cM_list_RA.txt", h=T, as.is=T)
rarthritis$Disease <- "Rheumatoid arthritis"
cad <- read.delim("WTCCC_DZ_allc_allp_GSeQTL_0.05cM_list_CAD.txt", h=T, as.is=T)
cad$Disease <- "Coronary artery disease"
crohns <- read.delim("WTCCC_DZ_allc_allp_GSeQTL_0.05cM_list_CD.txt", h=T, as.is=T)
crohns$Disease <- "Crohn's disease"
t1d <- read.delim("WTCCC_DZ_allc_allp_GSeQTL_0.05cM_list_T1D.txt", h=T, as.is=T)
t1d$Disease <- "Type 1 diabetes"
t2d <- read.delim("WTCCC_DZ_allc_allp_GSeQTL_0.05cM_list_T2D.txt", h=T, as.is=T)
t2d$Disease <- "Type 2 diabetes"

all <- rbind(hypertension, bipolar, rarthritis, cad, crohns, t1d, t2d)

# selection-enriched genes
genes <- read.table("../ld/hc_snps/genes_unique_hc.txt", h=F, as.is=T)
colnames(genes) <- "exGene"

find <- merge(genes, all, by="exGene")
find_pop <- merge(find, eqtl, by="exGene") #integrate with pop-specific data

#ceu <- filter(find_pop, Population=="CEU" & High_confidence=="Y")
#yri <- filter(find_pop, Population=="YRI" & High_confidence=="Y")
ceu <- filter(find_pop, Population=="CEU")
yri <- filter(find_pop, Population=="YRI")

ceu_unique <- ceu[!duplicated(ceu$eQTL), ]
yri_unique <- yri[!duplicated(yri$eQTL), ]

length(ceu_unique$eQTL)
[1] 1443
length(unique(ceu_unique$exGene))
[1] 32
length(yri_unique$eQTL)
[1] 725
> length(unique(yri_unique$exGene))
[1] 9

write.table(ceu_unique, "ceu_disease_eqtl.txt", col.names=T, row.names=F,
  quote=F, sep="\t")
write.table(yri_unique, "yri_disease_eqtl.txt", col.names=T, row.names=F,
  quote=F, sep="\t")
