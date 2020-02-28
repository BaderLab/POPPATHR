#' Sets up and runs GSEA using population-based FST
#' table containing GSEA results (Geneset, Size, ES, NES, NominalP, FDR, FWER)
#'
#' @param fst_file (char) path to file with SNP-level FST statistics
#' @param annotation_file (char) path to pathway definitions GMT file
#' @param snp2gene_file (char) path to file with snp2gene mappings. Output of
#' 		mapSNP2gene() (found in GWAS2Pathway).
#' @param SET_PERM (integer) set cycle of permutations to run
#' 		default=1000
#' @param SNP2GENE_DIST (integer) value for GSEA --distance.
#'		Max. distance between SNP and gene for their association in snp2gene_file
#' @param MIN_GENE (integer) value for GSEA --setmin.
#' 		Min. number of genes in a gene set to be considered
#' @param MAX_GENE (integer) value for GSEA --setmax.
#' 		Max. number of genes in a gene set to be considered
#' @param SET_SEED (integer) value for GSEA --seed
#'
#' @return none
#' @export
#'

setupGSEArun <- function(fst_file, annotation_file, snp2gene_file,
												 SET_PERM=SET_PERM, SNP2GENE_DIST=SNP2GENE_DIST,
												 MIN_GENE=MIN_GENE, MAX_GENE=MAX_GENE, SET_SEED=42,
												 output_folder) {
	# Setup and run GSEA
 	statOut <- sprintf("%s/gseaStatFile.txt", output_folder)
 	leOut	  <- sprintf("%s/gseaLEout.txt", output_folder)
 	resOut  <- sprintf("%s/results.txt", output_folder)

	# Run genotype-permuted GSEA
	cat("* Running genotype-based permutations\n\n")
	genoPermCom <- paste("calculate_gsea.pl %s %s --mapfile %s --cycle %i",
											 "--setstatfile %s --leout %s --setmin %i --setmax %i",
											 "--seed %i --distance %i --log %s/calculate_gsea.log")
	cmd <- sprintf(genoPermCom, fst_file, annotation_file, snp2gene_file,
								 SET_PERM, statOut, leOut, MIN_GENE, MAX_GENE,
								 SET_SEED, SNP2GENE_DIST, output_folder)
  system(cmd)

  # Combine GSEA results
	system(sprintf("combine_gsea.pl %s/calculate_gsea.log > %s/combined_res.log",
		  output_folder, output_folder))

  # Format results for output
	cat("* Formatting results for output...")
  dat <- read.delim(sprintf("%s/combined_res.log", output_folder),
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
	#write.xlsx(dat,file=sprintf("%s/results.xlsx", output_folder), col=TRUE, row=FALSE)

	# Write table of significant pathways to output directory
	write.table(dat[sigPaths1,], file=sprintf("%s/pathways_FDR0.1.txt", output_folder),
		sep="\t", col=TRUE, row=FALSE, quote=FALSE)
	write.table(dat[sigPaths2,], file=sprintf("%s/pathways_FDR0.05.txt", output_folder),
		sep="\t", col=TRUE, row=FALSE, quote=FALSE)
	cat(" done.\n")

  # Cleanup
  unlink(sprintf("%s/combined_res.log", output_folder))
}
