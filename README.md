# POPPATHR

POPPATHR is an R package that determines pathway-level SNP-SNP associations (coevolution) driven by population positive selection to better understand the evolution of human pathways.

## Prerequisites

R version > 3.6

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

### Installing

Install POPPATHR using devtools:

```
# Install devtools
install.packages("devtools")

# Install POPPATHR
devtools::install_github("BaderLab/POPPATHR")
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

Cytoscape needs to be downloaded from your browser using the link above.

### Data inputs

You will need three types of data to use POPPATHR:
- [SNP genotypes](#snp-genotypes)
- [Gene set annotations](#gene-set-annotations)
- [Functional genome annotations](#functional-genome-annotations)

#### SNP genotypes

In the paper, we used SNP genotype data from the International HapMap Project 3 (HM3). This data required a genome assembly conversion from hg18 to hg19, and needed to be converted to PLINK [binary format](http://www.cog-genomics.org/plink/1.9/formats) (bim, bed, fam).

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

We use a set of gene set annotations from the [Bader lab](http://download.baderlab.org/EM_Genesets/).

The set we use contains only gene sets from GO biological process excluding annotations that have evidence code IEA (inferred from electronic annotation), ND (no biological data available), and RCA (inferred from reviewed computational analysis) and all pathway resources.

```
# Enter POPPATHR annotations folder
# Assuming POPPATHR_DIR is defined as above
cd ${POPPATHR_DIR}/data/annotations

# Define path to file
home_folder="http://download.baderlab.org/EM_Genesets/"

# We used the file dated April 26 2016 in our paper
date="April_24_2016"
organism="Human"
key_type="symbol"
file_path="${organism}_GOBP_AllPathways_no_GO_iea_${date}_${key_type}.gmt"

# Download file
wget ${home_folder}/${date}/${organism}/${key_type}/${file_path}
```

#### Functional genome annotations

We assess the functionality of

## Basic usage

You've made it to this step, hurray!



## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags).

## Authors

* **Catherine Ross**

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
