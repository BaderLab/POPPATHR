# Prepare ChiP-peak bed file for FunciSNP biofeautre input

dataDir <- "/scratch/g/gbader/cmross/PopulationPathways/data"
# Dir with raw bed biofeature files
bioDir <- sprintf(paste0("%s/roadmap/chromhmmSegmentations/ChmmModels/",
                         "imputed12marks/jointModel/final/HSC_Bcell"), dataDir)

outDir <- sprintf("%s/for_funcisnp", bioDir)

#-------------------------------------
require(dplyr)

#===============================================================================
## Specific to GTEx eQTL bed files

#system(sprintf("gzip -d %s.gz", bioF))
#dat <- read.table(bioF, h=F, as.is=T)
#dat <- dat[,1:4] # get relevant columns
#dat <- dat[!is.na(as.numeric(as.character(dat[,1]))), ] # remove 'X' and 'Y' chr
#dat[,1] <- paste("chr", dat[,1], sep="") # add 'chr' to first column
#write.table(dat, file=sprintf("%s/FS_%s", outDir, basename(bioF)),
#            col=F, row=F, quote=F, sep="\t")

#===============================================================================
## Specific to Roadmap Epigenomics chromHmm peak bed files (mnemonics)
# Unzip gzipped bed files

#' @param bioDir (char) path to file containing raw bed files for FS input_snps
#' @param filter (logical) if TRUE filters files based on desired state
#' @param pattern (char) string indicating desired TS mnemonic
#' leave blank if filter == FALSE
#' eg. 'TssA'	= Active TSS, 'TxReg'	Transcribed & regulatory (Prom/Enh)
#' refer to http://egg2.wustl.edu/roadmap/web_portal/imputed.html#chr_imp

formatBio <- function(bioDir, filter=TRUE, pattern="TxReg") {

  # Decompress biofeature files
  system(sprintf("gzip -d %s/*", bioDir))
  bio <- list.files(path=bioDir, pattern="*.bed", full.names=T)
  cat(sprintf("Decompressed %i bed files in %s.\n", length(bio), bioDir))

  if (filter == TRUE) {
    # Read in all decompressed files
    bed_list <- lapply(bio, read.table)
    if (!file.exists(outDir)) dir.create(outDir)

    format_list <- list()
    for (i in 1:length(bed_list)) {
        cat(sprintf("\nFormatting %s biofeature file...", basename(bio[i])))
        format_list[[i]] <- filter(bed_list[[i]], grepl(pattern, V4, fixed=T))
        write.table(format_list[[i]],
                    file=sprintf("%s/%s_%s", outDir, pattern, basename(bio[i])),
                    col=F, row=F, quote=F, sep="\t")
        cat(sprintf(" file written to %s/%s_%s.\n", outDir, pattern,
                    basename(bio[i])))
    }
  }
  else {
    return(NULL)
  }

  num.org    <- sapply(bed_list, nrow)
  num.filter <- sapply(format_list, nrow)

  cat(sprintf("Pulled %i %s lines from %i total lines in %s peak file.\n",
              num.filter, pattern, num.org, basename(bio)))
}

#-------------------------------------
formatBio(bioDir, filter=TRUE, pattern="TxReg")
