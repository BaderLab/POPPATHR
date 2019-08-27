#' Sets up and runs GSEA using population-based FST
#' table containing GSEA results (Geneset, Size, ES, NES, NominalP, FDR, FWER)
#'
#' @param realF (char) path to file with real SNP association statistics
#' @param pathF (char) path to pathway definitions GMT file
#' @param snp2geneF (char) path to file with snp2gene mappings. Output of
#' 		mapSNP2gene() (found in GWAS2Pathway).
#' @param setPerm (integer) set cycle of permutations to run
#' 		default=1000
#' @param snp2genedist (integer) value for GSEA --distance.
#'		Max. distance between SNP and gene for their association in snp2geneF
#' @param minGene (integer) value for GSEA --setmin.
#' 		Min. number of genes in a gene set to be considered
#' @param maxGene (integer) value for GSEA --setmax.
#' 		Max. number of genes in a gene set to be considered
#' @param setSeed (integer) value for GSEA --seed
#'
#' @return none
#' @export
#'
setupGSEArun <- function(realF, pathF, snp2geneF,
												 setPerm=1000L, snp2genedist=500000L,
												 minGene=10L, maxGene=300L, setSeed=42L,
												 outDir) {
	# Setup and run GSEA
 	statOut <- sprintf("%s/gseaStatFile.txt", outDir)
 	leOut	  <- sprintf("%s/gseaLEout.txt", outDir)
 	resOut  <- sprintf("%s/results.txt", outDir)

	# Run genotype-permuted GSEA
	cat("* Running genotype-based permutations\n\n")
	genoPermCom <- paste("calculate_gsea.pl %s %s --mapfile %s --cycle %i",
											 "--setstatfile %s --leout %s --setmin %i --setmax %i",
											 "--seed %i --distance %i --log %s/calculate_gsea.log")
	cmd <- sprintf(genoPermCom, realF, pathF, snp2geneF,
								 setPerm, statOut, leOut, minGene, maxGene,
								 setSeed, snp2genedist, outDir)
  system(cmd)

  # Combine GSEA results
	system(sprintf("combine_gsea.pl %s/calculate_gsea.log > %s/combined_res.log",
		  outDir, outDir))

  # Format results for output
	cat("*Formatting results for output...")
  dat <- read.delim(sprintf("%s/combined_res.log", outDir),
		skip=1, h=FALSE, as.is=TRUE)
  dat[,1] <- sub("Geneset=","",	 dat[,1])
  dat[,2] <- sub("Size=","",		 dat[,2])
  dat[,3] <- sub("ES=","",		   dat[,3])
  dat[,4] <- sub("NES=","",		   dat[,4])
  dat[,5] <- sub("NominalP=","", dat[,5])
  dat[,6] <- sub("FDR=","",		   dat[,6])
  dat[,7] <- sub("FWER=","",		 dat[,7])
  colnames(dat) <- c("Geneset", "Size", "ES", "NES", "NominalP", "FDR", "FWER")
  dat <- dat[order(dat$FDR),]
  write.table(dat, file=resOut, sep="\t", col=TRUE, row=FALSE, quote=FALSE)
	cat(" done.\n")

	# Write output to excel format
	#write.xlsx(dat,file=sprintf("%s/results.xlsx", outDir), col=T, row=F)

	# Tally number of pathways at sig. thresholds
	sigPaths1 <- which(dat$FDR <= 0.1)
	sigPaths2 <- which(dat$FDR <= 0.05)
	cat("\n------SIGNIFICANT GSEA RESULTS------\n")
	cat("Number of significantly enriched pathways (FDR <= 0.1):",
			length(sigPaths1))
	cat("\nNumber of significantly enriched pathways (FDR <= 0.05):",
			length(sigPaths2))

	# Write table of significant pathways to output directory
	write.table(dat[sigPaths1,], file=sprintf("%s/pathways_FDR0.1.txt",
			outDir), sep="\t", col=TRUE, row=FALSE, quote=FALSE)
	write.table(dat[sigPaths2,], file=sprintf("%s/pathways_FDR0.05.txt",
			outDir), sep="\t", col=TRUE, row=FALSE, quote=FALSE)

  # Format results for Cytoscape enrichment map
	cat("\n\n*Formatting output for use in Cytoscape (Enrichment map)...")
  dat <- dat %>% mutate(Description=Geneset)
  dat <- subset(dat, select=c(Geneset, NominalP, FDR, Description))
  dat <- dat[c("Geneset", "Description", "NominalP", "FDR")]
  write.table(dat, file=sprintf("%s/results_forEM.txt", outDir),
        			sep="\t", col=TRUE, row=FALSE, quote=FALSE)
	cat(" done.\n")

  # Cleanup
  unlink(sprintf("%s/combined_res.log", outDir))
}
