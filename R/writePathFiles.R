#' Generates SNP lists per selection-enriched and unenriched pathway
#' as determined by GSEA
#'
#' @param genotypeFile (char) path to file with SNP genotype data (PLINK format).
#' @param resultsFile (char) path to files with GSEA results of two population
#' comparisons (i.e., CEU vs. YRI and CEU vs. LWK)
#' @param snp2geneFstFile (char) path to SNP-gene mapping with FST file.
#' @param famFile (char) path to PLINK population coded fam file.
#' @param emGroupFile (char) file path to write pathway groupings as determined
#' by EnrichmentMap.
#' @param enrichFDR (integer) FDR cutoff to select validated selection-enriched
#' pathways (default=0.05).
#' @param unenrichNES (integer) NES cutoff to select unenriched pathways
#' (default=1).
#' @param enrichFolder (char) path to directory to store output files
#' (PLINK files per selection-enriched gene set).
#' @param enrichEMFolder (char) path to directory to store output files
#' (PLINK files per selection-enriched pathway, grouped via AutoAnnotate).
#' @param unenrichFolder (char) path to directory to store output files
#' (PLINK files per unenriched gene set).
#' @export
writePathFiles <- function(genotypeFile, resultsFile,
	snp2geneFile, famFile, emGroupFile, enrichFDR=3, unenrichNES=1,
	enrichFolder, enrichEMFolder, unenrichFolder) {

	# Read in snp2geneFst file
	snpFst <- read.delim(snp2geneFstFile, h=TRUE, as.is=TRUE)

	# Read in GSEA results files formatted for EnrichmentMap
    cat(sprintf("* Reading in GSEA results file %s\n", resultsFile))
    resNames <- lapply(strsplit(resultsFile, "\\/"), "[", 2)
    resFiles <- lapply(resultsFile, function(x) read.delim(x, h=TRUE, as.is=TRUE))

    # Prepare data to combine
    for (i in seq_along(resFiles)) {
		resFiles[[i]] <- resFiles[[i]][,c("pathway", "pval", "padj", "NES", "leadingEdge")]
        colnames(resFiles[[i]])[2:ncol(resFiles[[i]])] <- paste(colnames(resFiles[[i]][2:ncol(resFiles[[i]])]), resNames[i], sep="_")
    }

	# Combine results from both population comparisons
    resComb <- purrr::reduce(resFiles, left_join, by="pathway")

	# Define function to grab pathway information (SNPs/genes)
	pathInfo <- function(pathwaySet, pathwayStat, outputFolder) {

		# Prepare output files
		pathwayOut <- sprintf("%s/results_%s.txt", outputFolder, pathwaySet)
		nameOut <- sprintf("%s/pathways_%s.txt", outputFolder, pathwaySet)
		statOut <- sprintf("%s/stats_%s.txt", outputFolder, pathwaySet)

		# Generate output directory for SNP/gene files per pathway set
		cat(sprintf("\n* Writing pathway files for %s pathways\n\n", pathwaySet))
		pathFolder <- sprintf("%s/pathwayFiles", outputFolder)
		if (!dir.exists(pathFolder)) dir.create(pathFolder)

		# Define pathways in specified pathway set
		cat(sprintf("* Getting table of %s pathways from GSEA results...", pathwaySet))
		if (pathwaySet == "enriched" | pathwaySet == "enrichedEM") {
			pathways <- filter(resComb, get(paste0("padj_", resNames[1])) <= 0.05)
			if (length(resFiles) > 1) {
				pathways <- filter(pathways, get(paste0("padj_", resNames[2])) <= 0.05)
			}
		}

		if (pathwaySet == "unenriched") {
			pathways <- filter(resComb, NES_CEU_YRI < 1)
			if (length(resFiles) > 1) {
				pathways <- filter(pathways, NES_CEU_LWK < 1)
			}
			# Subset unenriched set to be same size as enriched
			n_enrich <- length(readLines(file.path(enrichFolder, "pathways_enriched.txt")))
			pathways <- pathways[(nrow(pathways)-n_enrich+1):nrow(pathways),]
		}

		# Write out to file
		cat(sprintf(" writing out to %s\n", pathwayOut))
		write.table(pathways, file=pathwayOut, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

		# Define pathways to grab SNP/gene info for
		if (pathwaySet == "enrichedEM") {
			load(emGroupFile)
			pathwayNames <- names(emGroupList)
		} else {
			pathwayNames <- pathways[,"pathway"]
		}

		cat(sprintf("* Grabbing names for %i %s pathways...", length(pathwayNames), pathwaySet))
		cat(sprintf(" writing out to %s", nameOut))
		write.table(pathwayNames, file=nameOut, col=FALSE, row=FALSE, quote=FALSE)


		# Pull SNP/gene information per pathway
		statList <- list()
		for (i in seq_along(pathwayNames)) {
			print(i)

			if (pathwaySet == "enrichedEM") {
				# Leading edge genes for all pathways in EM group
				if (length(grep("\\+|\\(", emGroupList[[i]]))) {
					em_i <- grep(paste(emGroupList[[i]], collapse="|"), pathways[,"pathway"], fixed = TRUE)
				} else {
					em_i <- grep(paste(emGroupList[[i]], collapse="|"), pathways[,"pathway"])
				}
				pathway_i <- pathways[em_i, paste0("leadingEdge_", resNames[1])]
			} else {
				# Leading edge genes for single pathway
				pathway_i <- pathways[i, paste0("leadingEdge_", resNames[1])]
			}

			# Dataframe format
			pathway_i <- data.frame(gene=unlist(strsplit(pathway_i, " ")))
			pathway_i <- unique(pathway_i)

			# Join with snp/gene/fst info and get unique entries
			pathway_i <- left_join(pathway_i, snpFst, by="gene")
			pathway_i$pathway <- pathwayNames[i]

			# Store data in list
			statList[[i]] <- pathway_i

			# Replace spaces with underscore in pathway pathway_name_i string
			# And replace special characters with underscore for PLINK compatibility
			pathway_name_i <- gsub(" ", "_", pathwayNames[i], fixed=TRUE)
			pathway_name_i <- gsub('([[:punct:]])|\\s+', '_', pathway_name_i)

			# Write out separate lists for pathway SNPs/genes
			snpList <- as.data.frame(pathway_i[,"snp"])
			geneList <- as.data.frame(pathway_i[,"gene"])

			# Writing out SNP list
			snpFile <- file.path(sprintf("%s/%s.snps", pathFolder, pathway_name_i))
			cat(sprintf("\n* Generating list for SNPs in %s pathway...", pathway_name_i))
			write.table(snpList, file=snpFile, col=FALSE, row=FALSE, quote=FALSE)
			cat(" done.\n")

			# Writing out gene list
			geneFile <- file.path(sprintf("%s/%s.genes", pathFolder, pathway_name_i))
			cat(sprintf("* Generating list for genes in %s pathway...", pathway_name_i))
			write.table(geneList, file=geneFile, col=FALSE, row=FALSE, quote=FALSE)
			cat(" done.\n")

			# Subsetting PLINK files for pathway SNPs (for use in SNPassoc*.R)
			outFile <- gsub("\\..*", "", snpFile)
			str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s --extract %s", genotypeFile, genotypeFile, famFile, snpFile)
			str2 <- sprintf("--make-bed --allow-no-sex --out %s", outFile)
			cmd <- sprintf("%s %s", str1, str2)
			system(cmd)
		}

		# Write out table of all pathway stats
		statList <- bind_rows(statList)
		write.table(statList, file=statOut, col=TRUE, quote=FALSE, sep="\t")

		# Concatenate SNP/gene files into master file
		#concatFiles <- function(x) {
			# Define name for concatenated file
		#	cat(sprintf("* Combining all %s into master file...", x))
		#	outFile <- sprintf("%s/%s_%s", outputFolder, x, pathwaySet)

			# Grab all x files in pathway set, either SNPs or genes
		#	cmd <- sprintf("cat %s/*.%s > %s.txt", pathFolder, x, outFile)
		#	cat(sprintf(" writing out to %s.txt\n", outFile))
		#	system(cmd)

			# Read in file and filter for unique SNPs or genes
		#	cat(sprintf("* Subsetting for unique %s...", x))
		#	master <- read.table(sprintf("%s.txt", outFile), h=FALSE, as.is=TRUE)
		#	masterUnique <- as.data.frame(unique(master))
		#	cat(sprintf(" writing out to %s_unique.txt.\n", outFile))
		#	write.table(masterUnique, file=sprintf("%s_unique.txt", outFile), col=FALSE, row=FALSE, quote=FALSE)
		#}

		# Apply to generate master SNP and gene lists
		#invisible(mapply(concatFiles, c("snps", "genes")))
	}

	 # Run function for enriched and unenriched pathways
	 setList <- c("enriched", "enrichedEM", "unenriched")
	 nesList <- c(enrichFDR, enrichFDR, unenrichNES)
	 dirList <- c(enrichFolder, enrichEMFolder, unenrichFolder)
	 invisible(mapply(pathInfo, setList, nesList, dirList))
}
