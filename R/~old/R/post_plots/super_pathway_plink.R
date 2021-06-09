#workdir/[ld]/hc_em_groups
require(tools)

# hm3 data
bed <- "../../../../data/HM3_2010-05_phase3/HM3_pops_hg19.bed"
bim <- "../../../../data/HM3_2010-05_phase3/HM3_pops_hg19.bim"
fam <- "../../../../data/HM3_2010-05_phase3/HM3_pops_hg19_CEU_YRI.fam"

PLINK <- "../../../../../Software/plink_mac/plink"

snp_file <- list.files(pattern="*.snps")
outDir <- getwd()

for (i in 1:length(snp_file)) {
  str1 <- sprintf("%s --bed %s --bim %s --fam %s --extract %s",
                  PLINK, bed, bim, fam, snp_file[i])
  str2 <- sprintf("--make-bed --allow-no-sex --out %s/%s",
                  outDir, file_path_sans_ext(snp_file[i]))
  cmd <- sprintf("%s %s", str1, str2)
  system(cmd)
}
