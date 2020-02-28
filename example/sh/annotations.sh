#!/bin/bash
# Script to download annotation files used in manuscript
SCRIPT_DIR=$(pwd)
ANNO_DIR=$SCRIPT_DIR/anno
cd $ANNO_DIR

## Pathway annotation file
date="April_24_2016"
organism="Human"
key_type="symbol"
path_home="http://download.baderlab.org/EM_Genesets/"
path_file="${organism}_GOBP_AllPathways_no_GO_iea_${date}_${key_type}.gmt"
wget ${path_home}/${date}/${organism}/${key_type}/${path_file}

## BioGRID annotation file
biogridF="https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.4.163/BIOGRID-ORGANISM-3.4.163.tab2.zip"
wget ${biogridF}
unzip $( basename "$biogridF" )
find *BIOGRID* ! -name BIOGRID-ORGANISM-Homo_sapiens-3.4.163.tab2.txt -delete # keep only human interaction file

## Gene property annotations
### GWAS
gwasF="https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
gwasF_out="nhgri_gwas_hits.txt"
wget ${gwasF} -O ${gwasF_out}

### Drug-gene interactions
drugF="http://www.dgidb.org/data/interactions.tsv"
drugF_out="dridb_drug_gene_interactions.txt"
wget ${drugF} -O ${drugF_out}

### OMIM (disease-gene associations)
# NOTE: You will need a valid API key; request for access on the OMIM API page
omimF="https://data.omim.org/downloads/T-I-boz5TYq7zhXRSgU7-w/morbidmap.txt"
omimF_out="omim_disease_gene_interactions.txt"
wget ${omimF} -O ${omimF_out}

### GTex derived stats (expression variation across individuals)
# NOTE: large file after unpacking (5GB)
eqtlF="https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz"
eqtlF_out="gtex_eqtl_V7"
mkdir ${eqtlF_out}
wget ${eqtlF} -O ${eqtlF_out}.tar.gz
tar -zxvf ${eqtlF_out}.tar.gz -C ${eqtlF_out} --strip-components 1
gzip -d ${eqtlF_out}/*
rm ${eqtlF_out}.tar.gz
