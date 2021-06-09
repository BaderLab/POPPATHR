#' Sets up and runs GSEA using difference in minor allele frequency (MAF)
#' for use in runGSEA.R

#' @param realF (char) path to file with MAF diff SNP association statistics
#' @param pathFile (char) path to pathway definitions GMT file
#' @param snp2geneF (char) path to file with snp2gene mappings. Output of
#' mapSNP2gene() (found in GWAS2Pathway)
#' @param CALC_GSEA (char) path to calculate_gsea.pl
#' @param COMB_GSEA (char) path to combine_gsea.pl
#' @param statOut (char) path for GSEA --setstatfile
#' @param leOut (char) path for GSEA --leout
#' @param setSeed (integer) value for GSEA --seed
#' @return (char) table containing GSEA results
#' (Geneset, Size, ES, NES, NominalP, FDR, FWER)
#' @export

setupGSEArun <- function(realF, pathFile, snp2geneF,
												 snp2genedist=10000L, CALC_GSEA, COMB_GSEA,
												 statOut, leOut,
												 setSeed=42L, minGene=20L, maxGene=200L,
												 format4EM=TRUE) {
out <- tryCatch({

	# Outputs statement to log file saying if code is running on SciNet server
  # or on personal workstation
  x <- Sys.info()["nodename"]
  if (any(grep("^gpc", x))) cat("on scinet\n") else cat("on workstation\n")

  # Setup GSEA
  cat("Setting up GSEA...")
  str1 <- sprintf("perl %s %s %s", CALC_GSEA, realF, pathFile)
  str2 <- sprintf("--mapfile %s --setstatfile %s --leout %s",
  		 						snp2geneF, statOut, leOut)
	str3 <- sprintf("--setmin %i --setmax %i", minGene, maxGene)
  str4  <- sprintf("--seed %i --distance %i --log %s/calculate_gsea.log",
  		 						 setSeed, snp2genedist, resDir)
  cmd   <- sprintf("%s %s %s %s", str1, str2, str3, str4)
  system(cmd)
	cat(" done.\n")

  # Combine GSEA results
	cat("Combining GSEA results...")
	system(sprintf("perl %s %s/calculate_gsea.log > %s/combined_res.log",
		  					 COMB_GSEA, resDir, resDir))
	cat(" done.\n")

  # Format results for output
	cat("Formatting results for output...")
  dat <- read.delim(sprintf("%s/combined_res.log", resDir),skip=1,h=F,as.is=T)
  dat[,1] <- sub("Geneset=","",	dat[,1])
  dat[,2] <- sub("Size=","",		dat[,2])
  dat[,3] <- sub("ES=","",		dat[,3])
  dat[,4] <- sub("NES=","",		dat[,4])
  dat[,5] <- sub("NominalP=","",	dat[,5])
  dat[,6] <- sub("FDR=","",		dat[,6])
  dat[,7] <- sub("FWER=","",		dat[,7])
  colnames(dat) <- c("Geneset","Size","ES","NES","NominalP","FDR","FWER")
  dat <- dat[order(dat$FDR),]
  write.table(dat,file=sprintf("%s/results.txt", resDir),
        			sep="\t", col=T, row=F, quote=F)
	cat(" done.\n")

	# Write output to excel format
	write.xlsx(dat,file=sprintf("%s/results.xlsx", resDir), col=T, row=F)

	# Tally number of pathways at sig. threshold FDR < 0.05 and FDR < 0.1
	sigPaths1 <- which(dat$FDR < 0.1)
	sigPaths2 <- which(dat$FDR < 0.05)
	cat("\n------SIGNIFICANT GSEA RESULTS------\n")
	cat("Number of significantly enriched pathways with FDR < 0.1 =",
			length(sigPaths1))
	cat("\nNumber of significantly enriched pathways with FDR < 0.05 =",
			length(sigPaths2))

	# Write table of significant pathways to output directory
	write.table(dat[sigPaths1,], file=sprintf("%s/pathways_FDR0.1.txt",
							resDir), sep="\t", col=T, row=F, quote=F)
	write.table(dat[sigPaths2,], file=sprintf("%s/pathways_FDR0.05.txt",
							resDir), sep="\t", col=T, row=F, quote=F)

  # Format results for Cytoscape enrichment map
	if (format4EM == TRUE) {

		require(dplyr)
		cat("\n\nFormatting output for use in Cytoscape...")
	  dat <- dat %>% mutate( Description = Geneset ) # copy 'Geneset' column to 'Description'
	  dat <- subset(dat, select=c(Geneset, NominalP, FDR, Description))
	  dat <- dat[c("Geneset", "Description", "NominalP", "FDR")]
	  write.table(dat, file=sprintf("%s/results_forEM.txt", EMdir),
	        			sep="\t", col=T, row=F, quote=F)
		cat(" done.\n")

	}
  # Cleanup
  unlink(sprintf("%s/combined_res.log", resDir))

	})
  return(out)
  cat("...closing log.\n")
}
