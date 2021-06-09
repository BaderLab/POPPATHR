#!/bin/bash

# Run FunciSNP for all chunks of high-confidence GSEA pathway SNPs.
parentDir=/scratch/g/gbader/cmross/PopulationPathways
pathGroup=rand #rand, high_conf

# dir with SNP files
inDir=$parentDir/methods/3_funcisnp/infiles/$pathGroup/split_snps

# dir with bed files containing annotation
annoDir=$parentDir/data/roadmap/chromhmmSegmentations/ChmmModels/imputed12marks
annoDir=$annoDir/jointModel/final/HSC_Bcell/for_funcisnp

# output dir
outDir=$parentDir/methods/3_funcisnp/outfiles/roadmap_TxReg
outDir=$outDir/out_$(date +%F)_HSC_Bcell_TxReg_$pathGroup/all_tables

# where each job is stored
jobDir=$parentDir/methods/3_funcisnp/jobs

mkdir -p $outDir;

scriptdir=`pwd`
for f in $inDir/; do
    echo $f
    baseF=`basename $f`
    jobFile=${jobDir}/${baseF}.sh
    cat /dev/null > $jobFile;

    echo "#!/bin/bash" >> $jobFile
    echo "#PBS -l nodes=1:ppn=8,walltime=00:30:00" >> $jobFile;
    echo "#PBS -N $baseF" >> $jobFile;
    echo "module load intel/15.0.6 R/3.2.3" >> $jobFile   #updated R version
    echo "cd $scriptdir;"                   >> $jobFile;
    echo ""                                 >> $jobFile
    echo "Rscript FunciSNP_run2.R ${inDir}/${baseF} $annoDir $outDir" >> $jobFile

    chmod u+x $jobFile
    cd $jobDir
    qsub $jobFile
    cd $scriptdir
done
