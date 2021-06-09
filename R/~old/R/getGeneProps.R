# Integration of gene properties / features for GI queries
### Get gene properties / features
sink(sprintf("%s/enrichedGeneProperties.log", outDir))
cat(sprintf("Gene property directory: %s\n", outDir))

cat("\nGENE PROPERTY DOWNLOADS\n")
## GWAS associated genes
gwasF <- "https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
gwasF_out <- sprintf("%s/nhgri_gwas_hits.txt", outDir)

if (file.exists(gwasF_out)) {
  cat(sprintf("GWAS association file already exists! (downloaded: %s)\n", file.info(gwasF_out)$ctime))
} else {
  system(sprintf("wget %s -O %s", gwasF, gwasF_out))
  cat(sprintf("*Downloaded GWAS association file (%s) to %s\n", gwasF, basename(gwasF_out)))
}

## Targets of drugs (drug-gene interaction database)
drugF <- "http://www.dgidb.org/data/interactions.tsv"
drugF_out <- sprintf("%s/drug_gene_interactions.txt", outDir)

if (file.exists(drugF_out)) {
  cat(sprintf("Drug-gene interaction file already exists! (downloaded: %s)\n", file.info(drugF_out)$ctime))
} else {
  system(sprintf("wget %s -O %s", drugF, drugF_out))
  cat(sprintf("*Downloaded drug-gene interaction file (%s) to %s\n", drugF, basename(drugF_out)))
}

## Disease gene associations - OMIM
# NOTE: You will need a valid API key; request for access on the OMIM API page
omimF <- "https://data.omim.org/downloads/T-I-boz5TYq7zhXRSgU7-w/morbidmap.txt"
omimF_out <- sprintf("%s/omim_morbidmap.txt", outDir)

if (file.exists(omimF_out)) {
  cat(sprintf("Disease-gene association file exists already! (downloaded: %s)\n", file.info(omimF_out)$ctime))
} else {
  system(sprintf("wget %s -O %s", omimF, omimF_out))
  cat(sprintf("*Downloaded disease-gene association file (%s) to %s\n", omimF, basename(omimF_out)))
}

## GTex derived stats (variation in expression across individuals)
# NOTE: large file after unpacking (5GB)
eqtlF <- "https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz"
eqtlF_out <- sprintf("%s/gtex_eqtl_V7", outDir)

if (file.exists(eqtlF_out)) {
  cat(sprintf("GTEx-derived eQTL stat files exist already! (downloaded: %s)\n", file.info(eqtlF_out)$ctime))
} else {
  dir.create(eqtlF_out)
  system(sprintf("wget %s -O %s.tar.gz", eqtlF, eqtlF_out))
  system(sprintf("tar -zxvf %s.tar.gz -C %s --strip-components 1", eqtlF_out, eqtlF_out))
  system(sprintf("gzip -d %s/*", eqtlF_out))
  system(sprintf("rm %s.tar.gz", eqtlF_out))
  cat(sprintf("*Downloaded + unpacked GTEx derived stat file (%s) to %s.tar.gz\n", eqtlF, basename(eqtlF_out)))
}

### Check formatting for each file separately, then query each gene
cat("\nFILE STATS\n")
## GWAS hits (category-based gwas assocations from SP)
###EFO_0000618 nervous system disease
###EFO_0000589 metabolic disease
###EFO_0000319 cardiovascular
###EFO_0000540 immune system disease
###EFO_0005741 infectious disease
###EFO_0000311 cancer
###EFO_0005105 lipid or lipoprotein measure
gwasDir <- sprintf("%s/GWAShits_byCategory_190410_v4", outDir)
gwasF <- list.files(pattern="GWAShits.*txt", path=gwasDir, full.names=T)
gwas <- lapply(gwasF, function(x) read.delim(x, h=T, as.is=T))
name <- c("cancer", "cardiovascular", "immune system disease", "metabolic disease",
          "nervous system disease", "lipid or lipoprotein measure", "infectious disease")
gwas <- Map(cbind, gwas, disease_category=name) #Map -> mapply
gwas <- do.call("rbind", gwas)

gwas2 <- gwas[,c(16,1,15,3,9,13)]
# Filtering studies with sample size >= 1000 + p value <= 5x10-8
gwas2_filt <- filter(gwas2, sample_size >= 1000 & pvalue <= 5e-08)

gwas_df <- gwas2[,c(1:4)] # with rsID column
colnames(gwas_df) <- c("disease_category", "gene", "snp", "trait")
#gwas_df <- gwas2[,c(1,2,4)] # without rsID column
#colnames(gwas_df) <- c("disease_category", "gene", "trait")
gwas_df <- unique(gwas_df)

cat(sprintf(paste("GWAS\n",
                  "* %i unique gene-trait associations\n",
                  "* %i genes + %i traits\n"),
      nrow(gwas_df), length(unique(gwas_df[,2])), length(unique(gwas_df[,3]))))

## Drug-gene interactions
drug <- read.delim(drugF_out, h=T, as.is=T)
drug_df <- drug[,c(1,8)] # NOTE using 'drug_name' instead of 'drug_claim_primary_name'
drug_df <- unique(drug_df)
colnames(drug_df) <- c("gene", "dgidb_drug_name")
# remove rows with empty cells
drug_df$dgidb_drug_name <- gsub("^$|^ $", NA, drug_df$dgidb_drug_name)
drug_df <- na.omit(drug_df)
cat(sprintf(paste("DGIdb\n",
                  "* %i unique drug-gene interactions\n",
                  "* %i genes + %i drugs\n"),
      nrow(drug_df), length(unique(drug_df[,1])), length(unique(drug_df[,2]))))

## Disease gene associations
# NOTE OMIM file has empty lines that need to be removed
omim_new <- readLines(omimF_out)
omim_new <- omim_new[-c(1:3, 7629:length(omim_new))]
omimF_out <- substr(omimF_out, 0, nchar(omimF_out)-4)
write.table(omim_new,
            file=sprintf("%s_new.txt", omimF_out),
            col=F, row=F, quote=F, sep="\t")

omim <- read.delim(sprintf("%s_new.txt", omimF_out), h=T, as.is=T)
omim <- omim[,c(1,2)]

# Remove strange phenotype strings
toFind <- paste(c("\\{", "\\["), collapse="|")
omim <- omim[-grep(toFind, omim$X..Phenotype),]

# Get phenotype
omim_pheno <- strsplit(omim$X..Phenotype, ", ")
omim_pheno <- lapply(omim_pheno, '[', 1)
omim_pheno <- unlist(omim_pheno)

# Split delimited gene names column `Gene.Symbols` and insert as new rows (duplicate)
omim_split <- strsplit(omim$Gene.Symbols, ", ")
omim_df <- data.frame(omim_phenotype=rep(omim_pheno, sapply(omim_split, length)),
                      gene=unlist(omim_split))
omim_df <- unique(omim_df)
cat(sprintf(paste("OMIM\n",
                  "* %i unique disease-gene associations\n",
                  "* %i genes + %i diseases\n"),
      nrow(omim_df), length(unique(omim_df[,2])), length(unique(omim_df[,1]))))

## GTEx eQTLs
eqtlF <- list.files(pattern="*.egenes.txt$", path=eqtlF_out, full.names=T)
tissue <- substr(basename(eqtlF), 0, nchar(basename(eqtlF))-14)
eqtl <- lapply(eqtlF, function(x) read.delim(x, h=T, as.is=T))

# Add `tissue` column denoting eqtl origin
# NOTE that the *.egenes.txt.gz files contain data for all genes tested.
# To obtain the list of eGenes, select the rows with 'qval' ≤ 0.05.
for (i in seq_along(eqtl)) {
  eqtl[[i]]$gtex_eqtl_tissue <- tissue[i]
  eqtl[[i]] <- filter(eqtl[[i]], qval <= 0.05)
  eqtl[[i]] <- eqtl[[i]][,c(2,19,34)]
  colnames(eqtl[[i]])[1:2] <- c("gene", "snp")
}

eqtl_df <- do.call("rbind", eqtl)
eqtl_df <- unique(eqtl_df)
cat(sprintf(paste("GTEx\n",
                  "* %i unique cis-eQTL associations (qval <= 0.05)\n",
                  "* %i genes + %i tissues\n"),
      nrow(eqtl_df), length(unique(eqtl_df[,1])), length(unique(eqtl_df[,2]))))

### Read in query data
cat("\nSELECTION-ENRICHED DATA\n")
gene <- read.delim(geneF), h=F, as.is=T)
colnames(gene) <- "gene"

cat("*Getting pathway cluster info\n")
clust <- read.delim(clustF, h=T, as.is=T)
clust <- clust[,-c(4:5)]
colnames(clust)[4] <- "pathway_cluster"

gene_clust <- join(gene, clust)
cat(sprintf("*Integrating gene properties for %i genes...", length(unique(gene_clust$gene))))

### WORK BEGINS
geneProp_list <- list()

getGeneProperties <- function(prop_list) {
  for (i in seq_along(prop_list)) {
    cat(sprintf("Integrating %s property\n", prop_list[i]))
    prop_query <- join(gene_clust, get(prop_list[i]))
    prop_query <- unique(prop_query)
    prop_query <- na.omit(prop_query)
    prop_query[] <- lapply(prop_query, function(x) as.character(x))

    # Aggregate multiple rows based on common genes
    prop_query2 <- aggregate(.~gene, prop_query, function(x) toString(unique(x)))

    # Write to geneProp_list
    geneProp_list[[i]] <<- prop_query2
  }
}

outF <- sprintf("%s/selEnrich_geneProperties", outDir)

# Get gene properties
## Also run gwas without snps
gwas_df <- gwas_df[,-3]
prop_list <- c("gwas_df", "drug_df", "eqtl_df", "omim_df")
getGeneProperties(prop_list)

# Combine data for all queries and write to table
all_prop <- join_all(dfs=geneProp_list, by="gene", type="full")
all_prop <- all_prop[order(all_prop$disease_category),]
write.table(all_prop, sprintf("%s_pathwayCluster_GWAS_gene.txt", outF), col=T, row=F, quote=F, sep="\t")

cat(" done.\n")
sink()
