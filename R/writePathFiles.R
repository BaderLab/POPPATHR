#' Generates SNP lists per selection-enriched and unenriched pathway
#' as determined by GSEA
#'
#' @param genotype_file (char) path to file with SNP genotype data (PLINK format).
#' @param results_file (char) path to files with GSEA results.
#' 		Strutured to compare results of two population analyses
#' 		i.e., CEU vs. YRI and CEU vs. LWK.
#' @param gseaStat_file (char) path to GSEA statistics file.
#' @param snp2gene_file (char) path to SNP-gene mapping file.
#' @param fam_file (char) path to PLINK population coded fam file.
#' @param EM_group_file (char) file path to write pathway groupings as determined
#' 		by EnrichmentMap.
#' @param ENRICH_NES (integer) NES cutoff to select validated selection-enriched
#'		pathways (default=3)
#' @param UNENRICH_NES (integer) NES cutoff to select unenriched pathways
#'    (default=0.1)
#' @param enrich_folder (char) path to directory to store output files
#' 		(PLINK files per selection-enriched gene set)
#' @param enrichEM_folder (char) path to directory to store output files
#' 		(PLINK files per selection-enriched pathway, grouped via AutoAnnotate)
#' @param unenrich_folder (char) path to directory to store output files
#' 		(PLINK files per unenriched gene set)
#'
#' @return none
#' @export
#'

writePathFiles <- function(genotype_file, results_file, gseaStat_file,
													 snp2gene_file, fam_file, EM_group_file,
												 	 ENRICH_NES, UNENRICH_NES,
											   	 enrich_folder, enrichEM_folder, unenrich_folder) {

	# Merge GSEA results from both population comparisons and change column names
	# so the results from each comparison can be identified after merging
	gsea_results <- lapply(results_file, function(x) read.delim(x, h=TRUE, stringsAsFactors=FALSE))
	for (i in seq_along(gsea_results)) {
		colnames(gsea_results[[i]])[2:7] <- paste(colnames(gsea_results[[i]][, c(2:7)]), i, sep="_")
	}
	gsea_results_merge <- join_all(gsea_results, by="Geneset")

	# Define function to grab pathway information (SNPs/genes)
	pathInfo <- function(pathway_set, NES, output_folder) {
		# Generate output directory for SNP/gene files per pathway set
		path_folder <- sprintf("%s/pathway_files", output_folder)
		if (!dir.exists(path_folder)) { dir.create(path_folder) }

		# Define pathways in specified pathway set
		cat(sprintf("\n\n* Identifying %s pathways from GSEA results", pathway_set))
		if (pathway_set == "enriched" | pathway_set == "enrichedEM") {
			# Grab enriched pathways via defined NES cutoff
			cat(sprintf(" (enriched NES threshold = %g)\n", NES))
			pathways <- filter(gsea_results_merge, NES_1 >= NES & FDR_1 <= 0.05)
			if (length(gsea_results) > 1) { # if run, use second population comparison to validate results
				pathways <- filter(pathways, NES_2 >= NES & FDR_2 <= 0.05)
			}
		} else if (pathway_set == "unenriched") {
			# Grab unenriched pathways via defined NES cutoff
			cat(sprintf(" (unenriched NES threshold = %g)\n", NES))
			pathways <- filter(gsea_results_merge, NES_1 <= NES & NES_1 >= -NES)
			if (length(gsea_results) > 1) { # if run, use second population comparison to validate results
				pathways <- filter(pathways, NES_2 <= NES & NES_2 >= -NES)
			}
		}

		# Write out to file
		pathway_file <- sprintf("%s/results_%s.txt", output_folder, pathway_set)
		cat(sprintf("** Writing out list of %s pathways to %s\n", pathway_set, pathway_file))
		write.table(pathways, file=pathway_file, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

		# Define pathway names
		cat(sprintf("* Grabbing names for %i %s pathways\n", nrow(pathways), pathway_set))
		pathway_names <- pathways[,1] # pathway names
		name_file <- sprintf("%s/pathways_%s.txt", output_folder, pathway_set)
		cat(sprintf("** Writing out list of pathway names to %s\n", name_file))
		write.table(pathway_names, file=name_file, col=FALSE, row=FALSE, quote=FALSE)

		# Subset GSEA statistics by define pathway names
		## NOTE --colour-never flag specifically for OSX, issue with readLines
		## reading ANSI color codes from grep output
		cat("* Subsetting GSEA statistics file\n")
		stat_file <- sprintf("%s/gseaStatFile_%s.txt", output_folder, pathway_set)
		cmd <- sprintf("grep --colour=never -f %s %s > %s", name_file, gseaStat_file, stat_file)
		cat(sprintf("** Writing out file to %s\n", stat_file))
		system(cmd)

		# Get SNP/gene information for defined pathways
		cat("* Pulling pathway SNP/gene information\n")
		load(EM_group_file)
		no_col <- max(count.fields(stat_file)) # get max number of columns in file
		gsea_stat <- readLines(stat_file) # read in lines
		gsea_stat <- str_split_fixed(gsea_stat, "\t", no_col)
		gsea_stat <- t(gsea_stat) # transpose data

		# Define pathways to grab information from
		if (pathway_set == "enrichedEM") {
			pathway_grab <- EM_group_list
			pathway_name <- names(pathway_grab)
		} else {
			pathway_grab <- gsea_stat
			pathway_name <- pathway_grab[1,]
		}

		# Pull SNP/gene information per pathway
		for (i in seq_along(pathway_name)) {
			if (pathway_set == "enrichedEM") {
				# Remember, enrichedEM defines multiple pathway gene sets per pathway
				# So we are dealing with several gene sets in an EM-defined pathway in many cases
				pathway_i <- gsea_stat[,which(gsea_stat[1,] %in% pathway_grab[[i]])]
				# Remove pathway names from matrix
				if (!length(ncol(pathway_i))) {
					pathway_i <- pathway_i[-1]
				} else {
					pathway_i <- pathway_i[-1,]
				}
			} else {
				# Remove pathway name from matrix
				pathway_i <- pathway_grab[-1,i]
			}

			# Split original matrix into 3 columns (gene, snp, fst value)
			pathway_i <- as.data.frame(str_split_fixed(pathway_i, ",", 3))

			# Get unique SNPs/genes in pathway_i
			pathway_i[pathway_i == ""] <- NA
			pathway_i <- na.omit(pathway_i)
			pathway_i <- unique(pathway_i)

			# Replace spaces with underscore in pathway pathway_name_i string
			# And replace special characters with underscore for PLINK compatibility
			pathway_name_i <- gsub(" ", "_", pathway_name[i], fixed=TRUE)
			pathway_name_i <- gsub('([[:punct:]])|\\s+', '_', pathway_name_i)

			# Write out separate lists for pathway SNPs/genes
			snp_list <- as.data.frame(pathway_i[,2])
			gene_list <- as.data.frame(pathway_i[,1])

			# Writing out SNP list
			snp_file <- file.path(sprintf("%s/%s.snps", path_folder, pathway_name_i))
			cat(sprintf("\n* Generating list for SNPs in %s pathway...", pathway_name_i))
			write.table(snp_list, file=snp_file, col=FALSE, row=FALSE, quote=FALSE)
			cat(" done.\n")

			# Writing out gene list
			gene_file <- file.path(sprintf("%s/%s.genes", path_folder, pathway_name_i))
			cat(sprintf("* Generating list for genes in %s pathway...", pathway_name_i))
			write.table(gene_list, file=gene_file, col=FALSE, row=FALSE, quote=FALSE)
			cat(" done.\n")

			# Subsetting PLINK files for pathway SNPs (for use in SNPassoc*.R)
			out_file <- gsub("\\..*", "", snp_file)
			str1 <- sprintf("PLINK --bed %s.bed --bim %s.bim --fam %s --extract %s",
											genotype_file, genotype_file, fam_file, snp_file)
			str2 <- sprintf("--make-bed --allow-no-sex --out %s", out_file)
			cmd <- sprintf("%s %s", str1, str2)
			system(cmd)
		}

		# Concatenate SNP/gene files into master file
		concatFiles <- function(x) {
			cat(sprintf("* Combining all %s into master file\n", x))
			# Define name for concatenated file
			out_file <- sprintf("%s/%s_%s", output_folder, x, pathway_set)
			# Grab all x files in pathway set, either SNPs or genes
			cmd <- sprintf("cat %s/*.%s > %s.txt", path_folder, x, out_file)
			cat(sprintf("** Writing file to %s.txt\n", out_file))
			system(cmd)
			# Read in file and filter for unique SNPs or genes
			master <- read.table(sprintf("%s.txt", out_file), h=FALSE, as.is=TRUE)
			master_unique <- as.data.frame(unique(master))
			cat(sprintf("** Writing unique %s file to %s_unique.txt.\n", x, out_file))
			write.table(master_unique, file=sprintf("%s_unique.txt", out_file), col=FALSE, row=FALSE, quote=FALSE)
		}
		mapply(concatFiles, c("snps", "genes"))
	}

 # Run function for enriched and unenriched pathways
 set_list <- c("enriched", "enrichedEM", "unenriched")
 nes_list <- c(ENRICH_NES, ENRICH_NES, UNENRICH_NES)
 dir_list <- c(enrich_folder, enrichEM_folder, unenrich_folder)
 mapply(pathInfo, set_list, nes_list, dir_list)
}
