#!/bin/bash
# Script to download SNP genotype data (HM3)
SCRIPT_DIR=$(pwd)
SOFTWARE_DIR=$SCRIPT_DIR/software
GENO_DIR=$SCRIPT_DIR/genotypes
cd $GENO_DIR

# Path to HapMap 3 data plus chain file for hg18 to hg19 build conversion
genoF="HM3_2010_05_phase3"
home="https://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/2010-05_phaseIII"
chain="http://hgdownload.cse.ucsc.edu/goldenpath/hg18/liftOver/hg18ToHg19.over.chain"

# List of individuals belonging to each panel
popsF="relationships_w_pops_041510.txt"
wget ${home}/${popsF}

# Download SNP genotyping data + chain file
home=${home}/plink_format
mapF="hapmap3_r3_b36_fwd.consensus.qc.poly.map"
pedF="hapmap3_r3_b36_fwd.consensus.qc.poly.ped"

wget ${home}/${mapF}.gz ${home}/${pedF}.gz ${chain}.gz
gunzip *.gz

# Run liftOver
liftOverPlink.py -m ${mapF} -p ${pedF} -o ${genoF} \
                 -e ${SOFTWARE_DIR}/liftOver \
                 -c ${GENO_DIR}/$(basename -- ${chain})

# Convert to binary PLINK format
plink --file ${genoF} --make-bed --out ${genoF}
