#' Download 9 chromHMM files and create a file with the union of the intervals
#' via bedtools (for use in FunciSNP)
set -ex
cell_types=(Gm12878 H1hesc Hepg2 Hmec Hsmm Huvec K562 Nhek Nhlf)
dir=/media/catherine/DATAPART1/Data/PopulationPathways/data/chromHMM
outDir=${dir}/master

mkdir -p $dir; cd $dir;
mkdir -p $outDir;

for ct in "${cell_types[@]}"; do
    echo $ct;
    remote=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmm${ct}HMM.bed.gz
    F=$(basename $remote .gz)
    wget --quiet -O - $remote | zcat - | cut -f 1-4 | perl -pe 's/\d+_(.+)/$1/' > $F
done

# apt-get install bedtools (Ubuntu)
bedtools unionbedg -header -names "${cell_types[@]}" -i *.bed > ${outDir}/master.chromhmm.bed
