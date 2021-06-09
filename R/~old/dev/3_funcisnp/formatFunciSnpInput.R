#' Prepares input file for FunciSNP analysis
#' (https://github.com/Bioconductor-mirror/FunciSNP)

#' @param plinkF (char) path to file containing information for
#' chromosome, position, and rsID of studied variants (eg PLINK bim file)
#' @param snpDir (char) path to files with tag SNPs to filter main PLINK file
#' @return funciSnpF (char) file with 3 columns as per FunciSNP guidelines
#' 1-position (chrom:position), 2-rsID, 3-population (EUR, AFR, ASN or ALL)

#-----------------------------INPUT FILES---------------------------------------
dataDir   <- "/media/catherine/DATAPART1/Data/PopulationPathways"
pathGroup <- "rand" # high_conf vs rand

plinkF <- sprintf(paste0("%s/data/PNC/all/PNC_imputed_merged.CLEAN_FINAL_",
                         "5sd_noAxiom_CEU_ASW.bim"), dataDir)

snpDir <- sprintf("%s/input_snps/%s", dataDir, pathGroup)
outDir <- sprintf("%s/methods/3_funciSNP/infiles/%s", dataDir, pathGroup)

#--------------------------------WORK BEGINS------------------------------------
if (file.exists(outDir)) {
  cat("Directory exists! Not overwriting\n")
  Sys.sleep(3)
}

if (!file.exists(outDir)) dir.create(outDir)

#--------------------------------FUNCTION---------------------------------------
formatFunciSnpInput <- function(plinkF, snpDir,
                                pathString=pathGroup,
                                splitLines=TRUE, splitNum=20L) {

  # Read in original PLINK file with complete SNP genotype data
  all_snps <- read.table(plinkF, h=F, as.is=T)
  colnames(all_snps) <- c("chr", "snp", "bp", "pos")

  # Format input file per FunciSNP guidelines
  # First get the chromosome + bp position for each SNP in original data
  chr_pos <- as.data.frame(paste(all_snps$chr, all_snps$pos, sep=":"))
  chr_pos_rsID <- cbind(chr_pos, all_snps[2])

  # Filter SNPs based on lists of pathway SNPs
  path_snpsF <- list.files(path=snpDir, pattern="*.snps$", full.names=T)
  path_names <- gsub("\\%.*", "", basename(path_snpsF))

  # Generate input file for each pathway to analyze each separately
  for (i in 1:length(path_snpsF)) {
    cat(sprintf("\nFormatting %s SNPs for FunciSNP input...", path_names[i]))
    path_snps <- read.table(sprintf("%s", path_snpsF[i]), h=F, as.is=T)
    colnames(path_snps) <- "snp"
    df <- chr_pos_rsID[which(chr_pos_rsID$snp %in% path_snps$snp), ]
    ##copy df to add column with "EUR" and "AFR" to each SNP, eg
    ## chr1:pos1  rs1   EUR
    ## chr1:pos1  rs1   AFR
    df_copy <- df
    df$population <- "EUR"
    df_copy$population <- "AFR"
    FSinput <- rbind(df, df_copy)
    ##write files to output directory
    write.table(FSinput,
                file=sprintf("%s/%s.txt", outDir, path_names[i]),
                col=F, row=F, quote=F, sep="\t")
    cat(sprintf(" file written to %s dir.\n", outDir))
  }

  # Combine files to analyze all pathway SNPs together
  # Then split them to reduce computational load of FS analysis
  # (program bugs out if you send all SNPs in one file)
  if (splitLines == TRUE) {

    splitDir <- sprintf("%s/split_snps", outDir)
    allSnpF <- sprintf("%s/all_snps_unique.txt", snpDir) #file made via getSnpLists.R
                                                         #in bin/1_gsea dir

    if (!file.exists(splitDir)) dir.create(splitDir)

    dat <- read.table(allSnpF, h=F, as.is=T)
    colnames(dat) <- "snp"
    df <- chr_pos_rsID[which(chr_pos_rsID$snp %in% dat$snp), ]
    df_copy <- df
    df$population <- "EUR"
    df_copy$population <- "AFR"
    FSinput <- rbind(df, df_copy)

    write.table(FSinput,
		file=sprintf("%s/%s", splitDir, pathString),
		col=F, row=F, quote=F, sep="\t")

    cat("Splitting all input SNPs into smaller files...")
    str1 <- sprintf("split -l %i %s/%s", splitNum, splitDir, pathString)
    str2 <- sprintf("%s/%s_", splitDir, pathString)
    system(sprintf("%s %s", str1, str2))
    cat(sprintf(" files written to %s.\n", splitDir))

    # Remove original combined file
    unlink(sprintf("%s/%s", splitDir, pathString))
  }
}

#===============================================================================
formatFunciSnpInput(plinkF, snpDir, pathString=pathGroup,
                    splitLine=TRUE, splitNum=20L)
