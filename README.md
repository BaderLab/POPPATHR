# POPPATHR

POPPATHR is an R package that determines pathway-level SNP-SNP associations (coevolution) driven by population positive selection to better understand the evolution of human pathways.

## Prerequisites

R version > 3.5.0

Software:
- [GenGen](https://github.com/WGLab/GenGen.git)
- [PLINK v1.90](http://www.cog-genomics.org/plink2/)
- [Cytoscape](https://cytoscape.org/download.html)

Optional for genome assembly conversion:
- [liftOver](https://genome.sph.umich.edu/wiki/LiftOver)
- [liftOverPlink](https://github.com/sritchie73/liftOverPlink)

R packages:
- **CRAN**: tidyverse, dplyr, data.table, reshape2, gdata, RColorBrewer, cowplot, and argparse
- **Bioconductor**: GenomicRanges, snpStats, and RCy3

## Getting started

### Cloning

Clone the package locally from Github using git:

```
# Run on the command line (Terminal)
git clone https://github.com/BaderLab/POPPATHR.git
```

### Software

Downloading the software packages above may look something like this (on a Mac):

```
# Run the following lines on the command line (Terminal)

# Enter POPPATHR software folder
cd POPPATHR/
POPPATHR_DIR=$(pwd)
cd ${POPPATHR_DIR}/data/software

# Clone GenGen repository
git clone https://github.com/WGLab/GenGen.git

# Download and unzip PLINK
wget http://s3.amazonaws.com/plink1-assets/plink_mac_20200616.zip
unzip *.zip -d plink
rm *.zip

# Download liftOver binary and hg18 to hg19 chain file
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/liftOver ./
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/hgsql ./

# Clone liftOverPlink repository
git clone https://github.com/sritchie73/liftOverPlink.git

# Add executables to PATH
echo '' >> ~/.bash_profile
echo '# POPPATHR software executables' >> ~/.bash_profile
echo 'export PATH="'${POPPATHR_DIR}/data/software'/GenGen:$PATH"' >> ~/.bash_profile
echo 'export PATH="'${POPPATHR_DIR}/data/software'/plink:$PATH"' >> ~/.bash_profile
echo 'export PATH="'${POPPATHR_DIR}/data/software'/liftOverPlink:$PATH"' >> ~/.bash_profile

# Source file for changes to take effect
source ~/.bash_profile
```

Cytoscape needs to be downloaded from your browser using the link [above](#prerequisites).

### Data inputs

You will need two types of data to use POPPATHR:
- [SNP genotypes](#snp-genotypes)
- [Gene set annotations](#gene-set-annotations)

Plus optional data to dig deeper into your POPPATHR results:
- [Functional genome annotations](#functional-genome-annotations)

#### SNP genotypes

In the paper, we used population SNP genotype data from the International HapMap Project 3 (HM3). This data required a genome assembly conversion from hg18 to hg19, and needed to be converted to PLINK [binary format](http://www.cog-genomics.org/plink/1.9/formats) (bim, bed, fam).

```
# Enter POPPATHR genotypes folder
# Assuming POPPATHR_DIR is defined as above
cd ${POPPATHR_DIR}/data/genotypes

# Define paths to HM3 genotype data and chain file for hg18 to hg19 build conversion
genotype_file="HM3_2010_05_phase3"
home_folder="https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII"
chain_file="http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain"

# Download file defining individuals per population panel
population_file="relationships_w_pops_041510.txt"
wget ${home_folder}/${population_file}

# Re-define paths to HM3 genotype data
home_folder=${home_folder}/plink_format
map_file="hapmap3_r3_b36_fwd.consensus.qc.poly.map"
ped_file="hapmap3_r3_b36_fwd.consensus.qc.poly.ped"

# Download files and unpack
wget ${home_folder}/${map_file}.gz ${home_folder}/${ped_file}.gz ${chain_file}.gz
gunzip *.gz

# Run liftOverPlink to convert genome assembly
liftOverPlink.py -m ${map_file} -p ${ped_file} -o ${genotype_file} \
                 -e ${POPPATHR_DIR}/data/software/liftOver \
                 -c ${POPPATHR_DIR}/data/genotypes/$(basename -- ${chain_file})

# Convert to PLINK binary format
plink --file ${genotype_file} --make-bed --out ${genotype_file}
```

#### Gene set annotations

To identify selection-enriched pathways by gene set enrichment analysis (implemented via GenGen), we used a set of gene set annotations from the [Bader lab](http://download.baderlab.org/EM_Genesets/).

The particular set we used contains only gene sets from GO biological process excluding annotations that have evidence code IEA (inferred from electronic annotation), ND (no biological data available), and RCA (inferred from reviewed computational analysis) and all pathway resources.

```
# Enter POPPATHR annotations folder
# Assuming POPPATHR_DIR is defined as above
cd ${POPPATHR_DIR}/data/annotations

# Define path to file
home_folder="http://download.baderlab.org/EM_Genesets/"

# NOTE: we used the annotation file dated April 26 2016
date="April_24_2016"
organism="Human"
key_type="symbol"
file_path="${organism}_GOBP_AllPathways_no_GO_iea_${date}_${key_type}.gmt"

# Download file
wget ${home_folder}/${date}/${organism}/${key_type}/${file_path}
```

#### Functional genome annotations

To assess the functionality of the identified selection-enriched pathways, we integrated the pathway variants with various genomic annotation features: [GWAS traits and diseases](https://www.ebi.ac.uk/gwas/), [expression quantitative trait loci (eQTLs)](https://www.gtexportal.org/home/), [disease phenotypes](https://omim.org/), and [drug targets](http://www.dgidb.org/).

To download data disease phenotype data from OMIM, you will need to request access and obtain a valid API key (do so [here](https://www.omim.org/api)).

```
# Enter POPPATHR annotations folder
# Assuming POPPATHR_DIR is defined as above
cd ${POPPATHR_DIR}/data/annotations

# NHGRI-EBI GWAS file
gwas_file="https://www.ebi.ac.uk/gwas/api/search/downloads/alternative"
gwas_out="nhgri_gwas_hits.txt"
wget ${gwas_file} -O ${gwas_out}

# DGIdb drug-gene interaction file
drug_file="http://www.dgidb.org/data/interactions.tsv"
drug_out="dgidb_drug_gene_interactions.txt"
wget ${drug_file} -O ${drug_out}

# OMIM disease-gene association file
omim_file="https://data.omim.org/downloads/[YOUR-API-KEY]/morbidmap.txt"
omim_out="omim_disease_gene_interactions.txt"
wget ${omim_file} -O ${omim_out}

# GTEx genotype-tissue expression files
# NOTE: large file after unpacking (5GB)
eqtl_file="https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz"
eqtl_out="gtex_eqtl_V7"
mkdir ${eqtl_out}
wget ${eqtl_file} -O ${eqtl_out}.tar.gz
tar -zxvf ${eqtl_out}.tar.gz -C ${eqtl_out} --strip-components 1
gzip -d ${eqtl_out}/*
rm ${eqtl_out}.tar.gz
```

## Basic usage

You've made it to this step, hurray!

Running POPPATHR involves 3 parts that are split into the following scripts (used in this order):
- get_enrichment.R
- get_coevolution.R
- get_properties.R

The pipeline is currently set up to run as default on the inputs we used for our paper: SNP genotypes for two population comparisons (**CEU_YRI** and **CEU_LWK**) along with the annotation files outlined above. You can run the complete pipeline with these defaults by executing the following shell script on the command line: `sh run_POPPATHR.sh`

Otherwise, you can supply POPPATHR scripts with your own data. You can find detailed descriptions of all R script arguments in a handy command line interface:

```
Rscript get_enrichment.R --help

usage: get_enrichment.R [-h] [-v] [-p POPULATION_PAIR] [-g GENOTYPE_FILE]
                        [-t POPULATION_TABLE] [-a ANNOTATION_FILE]
                        [-r REFGENE_FILE] [-o OUTPUT_FOLDER]
                        [--SET_PERM SET_PERM] [--MIN_GENE MIN_GENE]
                        [--MAX_GENE MAX_GENE] [--SNP2GENE_DIST SNP2GENE_DIST]

POPPATHR: Population-based pathway analysis of SNP-SNP coevolution. This
script identifies selection-enriched pathways between two population cohorts.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Print verbose output
  -p POPULATION_PAIR, --population_pair POPULATION_PAIR
                        Names of two population cohorts to test (e.g., CEU_YRI
                        or CEU_LWK)
  -g GENOTYPE_FILE, --genotype_file GENOTYPE_FILE
                        Path to PLINK (bed, bim, fam formatted) SNP genotype
                        files [default data/genotypes/HM3_2010_05_phase3]
  -t POPULATION_TABLE, --population_table POPULATION_TABLE
                        Path to table defining population genotypes [default
                        data/genotypes/relationships_w_pops_041510.txt]
  -a ANNOTATION_FILE, --annotation_file ANNOTATION_FILE
                        Path to gmt file containing pathway annotations
                        [default data/annotations/Human_GOBP_AllPathways_no_GO
                        _iea_April_24_2016_symbol.gmt]
  -r REFGENE_FILE, --refgene_file REFGENE_FILE
                        Path to refGene genome annotation file [default
                        data/annotations/refGene.hg19.header.txt]
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Path to output folder [default output]
  --SET_PERM SET_PERM   Number of GSEA permutation cycles to run [default
                        10000]
  --MIN_GENE MIN_GENE   Minimum number of genes permitted in pathway gene set
                        [default 10]
  --MAX_GENE MAX_GENE   Maximum number of genes permitted in pathway gene set
                        [default 300]
  --SNP2GENE_DIST SNP2GENE_DIST
                        Maximum distance (bp) considered for SNP-to-gene
                        mapping [default 500000.0]
```

You will find two results folders in the **output** directory named by the population comparisons that were run. For example, if you ran `get_enrichment.R` on **CEU_YRI** and **CEU_LWK**, you will find two folders in **output** named accordingly with lots of data inside to peruse.

## Versioning

For the versions available, see the [tags on this repository](https://github.com/BaderLab/POPPATHR/tags).

## Authors

- **Catherine Ross**
- **Shraddha Pai**
