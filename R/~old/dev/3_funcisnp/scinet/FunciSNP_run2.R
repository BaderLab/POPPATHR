# Runs FunciSNP on a single input file, using annotation dir specified

require(FunciSNP)

# --------------------
# Input params

args    <- commandArgs(TRUE)

snpFile <- args[1] # /scratch/g/gbader/cmross/PopulationPathways/methods/1_gsea/PNC/CEU_ASW/out_161109_MAFdiff_CEU-ASW/post_analyses
annoDir <- args[2] # /scratch/g/gbader/cmross/PopulationPathways/data/cistrome
outDir  <- args[3] # /scratch/g/gbader/cmross/PopulationPathways/methods/3_funciSNP/outfiles

# directory where 1KGP vcf tabix files are located.
# the slash at the end of the path is required to work with FunciSNP.
KGP_dir <- "/scratch/g/gbader/cmross/PopulationPathways/anno/1000genomes/"

# --------------------
# Work begins

outFile <- sprintf("%s/%s_FStable.txt", outDir, basename(snpFile))
print(outFile)

fs  <- getFSNPs(snp.regions.file=snpFile,
                bio.features.loc=annoDir, par.threads=4L,
                built.in.biofeatures=FALSE,
                primary.server=KGP_dir)

if (nrow(fs@summary.data)<1) {
    cat("No overlapping features for any SNPS in list. Ending.\n")
} else {
    fs_anno <- FunciSNPAnnotateSummary(fs)
    write.table(fs_anno, file=outFile, sep="\t", col=T, row=F, quote=F)
}
