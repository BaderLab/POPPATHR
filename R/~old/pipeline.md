## PopulationPathways workflow

### 1) `recodeFAM.R`
* recodes PLINK fam file to case/control format by population
* supports running phenotype-based permutations (see `setupGSEArun.R`)

### 2) `popPCA.R`
* performs cluster analysis (PCA) between both studied populations via PLINK
--cluster
* outputs MDS plot (OPTIONAL)

### 3) `SNP2gene.R`
* integrates input SNP locations to their nearest gene
* code adapted from GWAS2pathway `map_snp2gene.R` (by Shraddha Pai)

### 4) `calcFST.R`
* calculates population-based differentiation in minor allele frequency (FST)
via PLINK --fst case-control
* outputs GSEA input SNP statistics file (SNP-FST)

### 5) `setupGSEArun.R`
* calculates pathway enrichment statistics via gene-set enrichment analysis
* selection-enriched pathway: enriched for high-FST genes in one population
compared to the other (i.e., enriched for positive selection)
* unenriched pathway: similar pathway-level FST between both populations
(i.e., not enriched for positive selection)
* code adapted from GWAS2pathway `GSEA_setup.R`

### 6) `getPathStats.R`
* generates SNP and gene lists per 'enriched' and 'unenriched' pathway as
determined by GSEA for downstream analysis

### 7) `LDstatsWPM.R`
* calculates extent and strength of coevolution *within* the pathways enriched
for positive selection
* uses inter-chromosomal linkage disequilibrium (LD) as a proxy for SNP-SNP
coevolution

### 7) `LDstatsBPM.R`
* calculates extent and strength of coevolution *between* pairs of pathways
enriched for positive selection
* uses inter-chromosomal linkage disequilibrium (LD) as a proxy for SNP-SNP
coevolution
