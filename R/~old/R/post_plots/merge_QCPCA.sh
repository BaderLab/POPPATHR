#!/bin/bash

#PBS -l nodes=1:ppn=8,walltime=00:30:00
#PBS -N merge_QPCA_info0.8

# QC of post-imputation merged, cleaned data.
# PCA against HM3 populations to ensure samples fall into major
# continental ancestry groups

module load intel R;

#cd /scratch/g/gbader/spai/BaderLab/PatientNetworks/PNC/impute

iodir=/scratch/g/gbader/spai/BaderLab/PNC/impute/merge_info0.8/all_chroms

PLINK=/home/g/gbader/spai/Other_Software/plink/plink
HM3=/scratch/g/gbader/spai/BaderLab/PNC/anno/hapmap3_pop/hg19/HM3_pops.hg19
HIGHLD=/scratch/g/gbader/spai/BaderLab/PNC/anno/high-LD-regions.txt
numSD=5 # num SD for ethnicity cutoff

### input plink file
inFile=${iodir}/PNC_imputed_merged.CLEAN
outPref=$inFile.forHM

###### 1. exclude snps in LD
###echo "* Pruning SNPs in LD"
### time $PLINK --bfile ${inFile} --exclude $HIGHLD \
###    --indep-pairwise 50 5 0.2 --make-bed --out ${outPref}_prune
######
####### 1b. compute cryptic relatedness, identify duplicates
######### Note prune.in contains SNPs retained after LD filtering
######### prune.out contains those excluded by LD filtering.
###$PLINK --bfile ${outPref}_prune \
###    --extract ${outPref}_prune.prune.in \
###    --genome --out ${outPref}_IBS
###Rscript plink_IBS.R ${outPref}_IBS.genome
###
####
####### 2. extract snps that are in linkage equilibrium
####### in prep for merging
###echo "  * Extract SNPs in linkage equilibrium from our data"
###time $PLINK --bfile $inFile --extract ${outPref}_prune.prune.in --make-bed --out ${outPref}_DATA_pruned
####
####
####### 3. now extract same SNPs from the HM3 reference data
#####echo "  * Extract our SNPs from HM3 data"
###$PLINK --bfile $HM3 --extract ${outPref}_prune.prune.in --filter-founders --make-bed --out ${iodir}/HapMapMini
###
###exit 0
##### 4. merge with HM3
##### This merge happens twice. First we identify snps that
##### don't align
##### and that we attempt to flip. We also identify and
##### exclude symmetric SNPs
##### Once the failures are resolved we do a "proper" second
##### merge.
###
##### 4a. first merge
###$PLINK --bfile ${outPref}_DATA_pruned --bmerge ${iodir}/HapMapMini.bed ${iodir}/HapMapMini.bim ${iodir}/HapMapMini.fam --geno 0.03 --make-bed --out ${iodir}/4Flips
###### 4b. compile symmetric snps
###echo "  * Compiling symmetric snps"
###symSNP=${iodir}/syms.txt
###cat /dev/null > $symSNP;
###awk '($5=="C" && $6 == "G"){print $2}' ${inFile}.bim >> $symSNP;
###awk '($5=="G" && $6 == "G"){print $2}' ${inFile}.bim >> $symSNP;
###awk '($5=="A" && $6 == "T"){print $2}' ${inFile}.bim >> $symSNP;
###awk '($5=="T" && $6 == "A"){print $2}' ${inFile}.bim >> $symSNP;
###
###
####### 4c. second proper merge
###cat $symSNP ${iodir}/HMDATA-merge.missnp > ${iodir}/problem_snps.txt
###echo "  * Merge using flip for missing snps"
###$PLINK --bfile  ${outPref}_DATA_pruned \
###    --flip ${iodir}/4Flips-merge.missnp \
###    --exclude ${iodir}/problem_snps.txt --make-bed --out ${outPref}_DATA_Align
###
###echo "	* Merge with HapMapv3"
###$PLINK --bfile ${outPref}_DATA_Align \
###    --bmerge ${iodir}/HapMapMini.bed ${iodir}/HapMapMini.bim ${iodir}/HapMapMini.fam \
###    --geno 0.03 --make-bed --out ${iodir}/HMDATA
###exit 0

###### if HMDATA-merge.missnp is detected
###$PLINK --bfile ${outPref}_DATA_Align --exclude ${iodir}/HMDATA-merge.missnp --make-bed \
###    --out ${outPref}_DATA_Align2
###$PLINK --bfile ${outPref}_DATA_Align2 \
### --bmerge ${iodir}/HapMapMini.bed ${iodir}/HapMapMini.bim ${iodir}/HapMapMini.fam \
###    --geno 0.03 --make-bed --out ${iodir}/HMDATA
###
###exit 0
##### 5. Now compute pairwise distances with HM3 samples
###echo "	* Generate pairwise IBS with HapMapv3"
###$PLINK --bfile ${iodir}/HMDATA --genome --out ${iodir}/IBS_HM_DATA
###
####### 6. Cluster with HM3 and plot.
###$PLINK --bfile ${iodir}/HMDATA \
###    --read-genome ${iodir}/IBS_HM_DATA.genome \
###    --cluster --mds-plot 3 --out ${iodir}/MDS

echo "If you see that patients mostly separate with HM3 groups all is well!"
Rscript plink_plotMDS.R ${outPref}_DATA_pruned.fam ${iodir}/MDS.mds $outPref $numSD
#Rscript plink_plotMDS.R ${outPref}_DATA_pruned.fam ${iodir}/MDS.mds $outPref   3
