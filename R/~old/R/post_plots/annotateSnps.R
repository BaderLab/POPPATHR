#' Annotate functionality of SNPs using AnnoVar
#' (http://annovar.openbioinformatics.org/en/latest/#annovar-documentation)

#' @param dbsnpFile (char) path to AnnoVar SNP annotation file
#' @param TABLE_AV (char) path to table_annovar.pl
#' @param ANNO_VAR (char) path to annotate_variation.pl
#' @param annodb (char) path to AnnoVar annotation directory
#' (eg. for purposes of this function, using hg19 build annotation)
#' @param annoBuild (char) string indicating genome build
#' @param dbString (char) string indicating the databases that will be used
#' to annotate each SNP provided in snpList
#' @param opString (char) string indicating which operation to use for each
#' database in dbString (eg. g=gene-based, r=region-based, f=filter-based)
#' @return (char) files containing AnnoVar annotations
#' @export

annotateSnps <- function(hcInDir, lcInDir, dbsnpFile, annodb,
                         annoBuild="hg19", TABLE_AV, ANNO_VAR,
                         dbString=paste0("refGene,cytoBand,genomicSuperDups,",
                                         "avsnp144,ljb26_all,1000g2015aug_all"),
                         opString="g,r,r,f,f,f") {

  doStuff <- function(snpDir, outDir) {

  snps <- list.files(path=snpDir, pattern="*.snps", full.names=T)

  for (i in 1:length(snps)) {
    message(sprintf("Creating AnnoVar input file for %s SNPs...", snps[i]))
    system(sprintf("fgrep -w -f %s %s > %s/%s_snplist.avinput",
                   snps[i], dbsnpFile, outDir, snps[i]))
    message(" done.\n")

    # Annotate variants in the avinput file
    message("AnnoVar processes beginnning...")
    str1 <- sprintf("perl %s %s/%s_snplist.avinput %s",
                    TABLE_AV, outDir, snps[i], annodb)
    str2 <- sprintf("-buildver %s -out %s/%s_snp_anno",
                    annoBuild, outDir, snps[i])
    str3 <- sprintf("-remove -protocol %s -operation %s", dbString, opString)
    str4 <- sprintf("-nastring .")
    cmd  <- sprintf("%s %s %s %s", str1, str2, str3, str4)
    system(cmd)

    # Run gene-based annotations (ie. classify them as intergenic, intronic,
    # non-synonymous SNP, frameshift deletion, large-scale duplication, etc.)
    str1 <- sprintf("perl %s -geneanno", ANNO_VAR)
    str2 <- sprintf("-buildver %s %s/%s_snplist.avinput %s",
                    annoBuild, outDir, snps[i], annodb)
    cmd <- sprintf("%s %s", str1, str2)
    system(cmd)

    # Run region-based annotations (ie. annotate variants that fall within
    # conserved regions)
    str1 <- sprintf("perl %s -regionanno", ANNO_VAR)
    str2 <- sprintf("-build %s -out %s/%s_phastCons100way",
                    annoBuild, outDir, snps[i])
    str3 <- sprintf("-dbtype phastConsElements100way %s/%s_snplist.avinput %s",
                    outDir, snps[i], annodb)
    cmd <- sprintf("%s %s %s", str1, str2, str3)
    system(cmd)

    # Run filter-based annotations (ie. identify a subset of variants input file
    # that are not observed in 1000G version 2014 Oct and those that are observed
    # with allele frequencies)
    #cat("Running filter-based annotations...")
    #str1 <- sprintf("perl %s -filter", ANNO_VAR)
    #str2 <- sprintf("-dbtype 1000g2014oct_all -buildver %s %s/snplist.avinput %s",
    #                annoBuild, outDir, annodb)
    #cmd <- sprintf("%s %s", str1, str2)
    #system(cmd)
    #cat(" done.\n")

    # Print summary of SNP annotations
    cat(sprintf("\n------Annotation summary for %s SNPs------\n", snps[i]))
    dat <- read.delim(sprintf("%s/%s_snp_anno.hg19_multianno.txt",
                      outDir, snps[i]), h=T, as.is=T)

    # Write output to excel format
    write.xlsx(dat, file=sprintf("%s/%s_snp_anno.hg19_multianno.xlsx",
               outDir, snps[i]), col=T, row=F)

     # Pie chart of SNP functional annotations (all annotated SNPs)
     annoTable <- table(dat$Func.refGene)
     annoTabledf <- as.data.frame(annoTable)
     pred_anno <- annoTabledf$Var1
     freq_anno <- annoTabledf$Freq
       freq_anno_pct <- round(freq_anno/sum(freq_anno) * 100)
       freq_anno_pct <- paste(freq_anno_pct, "%", sep="")
       freq_anno_pct <- sprintf("(%s)", freq_anno_pct)
     plot_anno_lbls <- paste(pred_anno, freq_anno_pct, sep=" ")

     png(file=sprintf("%s/%s_function_piechart.png",
         outDir, snps[i]))
     pie(annoTable, labels=plot_anno_lbls,
         col=rainbow_hcl(length(annoTable)),
         main=paste0(sprintf("Proportion of %s variant function", snps[i])))
     dev.off()

    # Rearrange data to get columns in the following order: "avsnp144 (rsID)",
    # "Func.refGene", "ExonicFunc.refGene", "SIFT_score", "SIFT_pred",
    # "Polyphen2_HDIV_score", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_score",
    # "Polyphen2_HVAR_pred", etc..
    dat2 <- dat[c(13,6,9,14:ncol(dat))]
    cat(sprintf("Total number of annotated SNPs: %i\n", nrow(dat2)))
    cat("Function of all annotated SNPs in refGene database:")
    print(table(dat2[2]))

    # Determining predictions
    cat("Functional predictions per annotation database:\n")
    # PolyPhen2 HDIV prediction
    type_of_change <- data.frame()    # Initialize empty data frame
    exonic <- dat2[which(dat2$Func.refGene == "exonic"),]

    polyphenPred <- exonic$Polyphen2_HDIV_pred
    for (i in 1:nrow(exonic)) {
      if(polyphenPred[i] == ".") {
        type_of_change[i, "Polyphen2_HDIV_pred"] <- "missing"
      }
      else if(polyphenPred[i] %in% "D") {
        type_of_change[i, "Polyphen2_HDIV_pred"] <- "probably_damaging"
      }
      else if(polyphenPred[i] %in% "P") {
        type_of_change[i, "Polyphen2_HDIV_pred"] <- "possibly_damaging"
      }
      else if(polyphenPred[i] %in% "B") {
        type_of_change[i, "Polyphen2_HDIV_pred"] <- "benign"
      }
      else {
        type_of_change[i, "Polyphen2_HDIV_pred"] <- "other"
      }
    }

    # PolyPhen2 HVAR prediction
    polyphenPred2 <- exonic$Polyphen2_HVAR_pred
    for (i in 1:nrow(exonic)) {
      if(polyphenPred2[i] == ".") {
        type_of_change[i, "Polyphen2_HVAR_pred"] <- "missing"
      }
      else if(polyphenPred2[i] %in% "D") {
        type_of_change[i, "Polyphen2_HVAR_pred"] <- "probably_damaging"
      }
      else if(polyphenPred2[i] %in% "P") {
        type_of_change[i, "Polyphen2_HVAR_pred"] <- "possibly_damaging"
      }
      else if(polyphenPred2[i] %in% "B") {
        type_of_change[i, "Polyphen2_HVAR_pred"] <- "benign"
      }
      else {
        type_of_change[i, "Polyphen2_HVAR_pred"] <- "other"
      }
    }

    # SIFT prediction
    SIFTpred <- exonic$SIFT_pred
    for (i in 1:nrow(exonic)) {
     if(SIFTpred[i] == ".") {
       type_of_change[i, "SIFT_pred"] <- "missing"
     }
     else if(SIFTpred[i] == "D") {
       type_of_change[i, "SIFT_pred"] <- "damaging"
     }
     else if(SIFTpred[i] == "T") {
       type_of_change[i, "SIFT_pred"] <- "tolerated"
     }
     else {
       type_of_change[i, "SIFT_pred"] <- "other"
     }
    }

    # LRT prediction
    LRTpred <- exonic$LRT_pred
    for (i in 1:nrow(exonic)) {
      if(LRTpred[i] == ".") {
        type_of_change[i, "LRT_pred"] <- "missing"
      }
      else if(LRTpred[i] == "D") {
        type_of_change[i, "LRT_pred"] <- "deleterious"
      }
      else if(LRTpred[i] == "N") {
        type_of_change[i, "LRT_pred"] <- "neutral"
      }
      else if(LRTpred[i] == "U") {
        type_of_change[i, "LRT_pred"] <- "unknown"
      }
    }

     # MutationTaster prediction
     MTpred <- exonic$MutationTaster_pred
     for (i in 1:nrow(exonic)) {
       if(MTpred[i] == ".") {
         type_of_change[i, "MutationTaster_pred"] <- "missing"
       }
       else if(MTpred[i] == "A") {
         type_of_change[i, "MutationTaster_pred"] <- "disease_causing_automatic"
       }
       else if(MTpred[i] == "D") {
         type_of_change[i, "MutationTaster_pred"] <- "disease_causing"
       }
       else if(MTpred[i] == "N") {
         type_of_change[i, "MutationTaster_pred"] <- "polymorphism"
       }
       else if(MTpred[i] == "P") {
         type_of_change[i, "MutationTaster_pred"] <- "polymorphism_automatic"
       }
       else {
         type_of_change[i, "MutationTaster_pred"] <- "other"
       }
     }

     # MutationAssessor prediction
     MApred <- exonic$MutationAssessor_pred
     for (i in 1:nrow(exonic)) {
       if(MApred[i] == ".") {
         type_of_change[i, "MutationAssessor_pred"] <- "missing"
       }
       else if(MApred[i] == "H") {
         type_of_change[i, "MutationAssessor_pred"] <- "high_functional"
       }
       else if(MApred[i] == "M") {
         type_of_change[i, "MutationAssessor_pred"] <- "medium_functional"
       }
       else if(MApred[i] == "L") {
         type_of_change[i, "MutationAssessor_pred"] <- "low_non-functional"
       }
       else if(MApred[i] == "N") {
         type_of_change[i, "MutationAssessor_pred"] <- "neutral_non-functional"
       }
       else {
         type_of_change[i, "MutationAssessor_pred"] <- "other"
       }
     }

     # FATHMM prediction
     FATHMMpred <- exonic$FATHMM_pred
     for (i in 1:nrow(exonic)) {
       if(FATHMMpred[i] == ".") {
         type_of_change[i, "FATHMM_pred"] <- "missing"
        }
       else if(FATHMMpred[i] == "D") {
          type_of_change[i, "FATHMM_pred"] <- "damaging"
        }
       else if(FATHMMpred[i] == "T") {
         type_of_change[i, "FATHMM_pred"] <- "tolerated"
       }
       else {
         type_of_change[i, "FATHMM_pred"] <- "other"
       }
     }

     # GERP++_RS prediction
     GERPpred <- as.numeric(exonic$GERP)
     for (i in 1:nrow(exonic)) {
       if(GERPpred[i] %in% NA) {
         type_of_change[i, "GERP"] <- "missing"
        }
       else if(GERPpred[i] > 0) {
          type_of_change[i, "GERP"] <- "conserved"
         }
       else if(GERPpred[i] < 0) {
          type_of_change[i, "GERP"] <- "not_conserved"
         }
       else {
          type_of_change[i, "GERP"] <- "other"
         }
      }

    # Print table of AnnoVar predictions for each exonic SNV to log
    exonic_SNVfunc <- exonic[c(1,3)]
    colnames(exonic_SNVfunc) <- c("Exonic_SNV", "Function")
    exonic_predTable <- cbind(exonic_SNVfunc, type_of_change)
    print(exonic_predTable)
  }

sink(sprintf("%s/annotateSnps.log", outDir))

  }

  doStuff(snpDir=hcInDir, outDir=)
  doStuff(snpDir=lcInDir, outDir=)

}
