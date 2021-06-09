# [workdir]/ld/nedelec_qtls

dat <- read.delim("putative_selected_qtl_nedelec_2016.txt", h=T, as.is=T)

sel_genes <- read.table("../genes_unique_hc.txt", h=F, as.is=T)
colnames(sel_genes) <- "HUGO_gene_id"
sel_snps <- read.table("../snps_unique_hc.txt", h=F, as.is=T)
colnames(sel_snps) <- "SNP.ID"

gene_overlap <- merge(dat, sel_genes, by="HUGO_gene_id")
nrow(gene_overlap) #36
length(unique(gene_overlap$HUGO_gene_id)) #32

snp_overlap <- merge(dat, sel_snps, by="SNP.ID")
nrow(snp_overlap) #0
