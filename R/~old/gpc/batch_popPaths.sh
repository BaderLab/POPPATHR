#!/bin/bash
# Runs PopulationPathways pipeline on Scinet GPC
parentDir=/scratch/g/gbader/cmross/PopulationPathways

# User parameters
data=HM3           #HM3
pop1=CEU           #CEU ASW
pop2=ASW           #LWK YRI
geno=HM3_pops_hg19
#geno=1kg_phase1_all_autosome_clean_ldpruned
pops=pop_info.txt        #for HM3
#pops=igsr_samples.tsv   #for 1KG

# Output directories
outDir=$parentDir/res/out_$(date +%F)_${geno}_${pop1}-${pop2}_MAF_20-200gene_sans_HLA
plinkDir=$outDir/freq
gseaDir=$outDir/gsea
resDir=$gseaDir/pathway_analysis

# Input files
genoF=$parentDir/data/$data/$geno
popInfo=`dirname ${genoF}`/${pops}
realFAM=${genoF}_${pop1}_${pop2}

# Software and annotation files
PLINK=/home/g/gbader/cmross/software/plink_linux_1.90/plink
CALC_GSEA=/home/g/gbader/cmross/software/GenGen-1.0.1/calculate_gsea.pl
COMB_GSEA=/home/g/gbader/cmross/software/GenGen-1.0.1/combine_gsea.pl

geneFile=$parentDir/anno/refGene/refGene.hg19.header.txt
#pathFile=$parentDir/anno/baderlab/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol.gmt
pathFile=$parentDir/anno/baderlab/Human_GOBP_AllPathways_no_GO_iea_April_24_2016_symbol_sans_HLA.gmt

# Where each job is stored
jobDir=$parentDir/res/jobs

# Work begins
mkdir -p $outDir;
mkdir -p $gseaDir;
mkdir -p $plinkDir;
mkdir -p $resDir;

scriptdir=`pwd`

baseF=`basename $realFAM`
jobFile=${jobDir}/${baseF}.sh
cat /dev/null > $jobFile;

echo "#!/bin/bash" >> $jobFile
echo "#PBS -l nodes=1:ppn=8,walltime=01:30:00" >> $jobFile; #add :m32g to request 32G node
echo "#PBS -N $baseF"                   >> $jobFile;
echo "module load intel/15.0.6 R/3.2.3" >> $jobFile   #updated R version
echo "cd $scriptdir;"                   >> $jobFile;
echo ""                                 >> $jobFile
echo "Rscript recodeFAM.R $data $pop1 $pop2 $genoF $popInfo \
                          `basename $realFAM` $outDir"          >> $jobFile
echo ""                                                         >> $jobFile
#echo "Rscript popPCA.R $genoF $realFAM $PLINK $pop1 $pop2 \
#                       ${pcaDir}/${pop1}_${pop2}"               >> $jobFile
#echo ""                                                         >> $jobFile
echo "Rscript SNP2gene.R ${genoF}.bim $geneFile \
                         ${gseaDir}/snp2gene.txt"               >> $jobFile
echo ""                                                         >> $jobFile
echo "Rscript calcMAFdiff.R $genoF ${realFAM}.fam $PLINK $plinkDir"   >> $jobFile
echo ""
echo "Rscript setupGSEArun.R ${plinkDir}/markerMAF.txt $pathFile ${gseaDir}/snp2gene.txt \
                             $CALC_GSEA $COMB_GSEA $resDir"     >> $jobFile

chmod u+x $jobFile
cd $jobDir
qsub $jobFile -q debug
cd $scriptdir
