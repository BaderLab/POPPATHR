#!/bin/bash
# Script to download bioinformatics software required for PopPaths pipeline
SCRIPT_DIR=$(pwd)
SOFTWARE_DIR=$SCRIPT_DIR/software
cd $SOFTWARE_DIR

# liftOver
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/liftOver ./
# hgsql
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/macOSX.x86_64/hgsql ./
# liftOverPLINK
git clone https://github.com/sritchie73/liftOverPlink.git
# PLINK
wget http://s3.amazonaws.com/plink1-assets/plink_mac_20190617.zip
unzip *.zip -d plink
rm *.zip
# GenGen
git clone https://github.com/WGLab/GenGen.git

# Add executables to PATH
echo '' >> ~/.bash_profile
echo '# Setting PATH for PopPaths software executables' >> ~/.bash_profile
echo 'export PATH="'${SOFTWARE_DIR}'/liftOverPlink:$PATH"' >> ~/.bash_profile
echo 'export PATH="'${SOFTWARE_DIR}'/plink:$PATH"' >> ~/.bash_profile
echo 'export PATH="'${SOFTWARE_DIR}'/GenGen:$PATH"' >> ~/.bash_profile
