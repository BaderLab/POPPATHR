% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/writePathFiles.R
\name{writePathFiles}
\alias{writePathFiles}
\title{Generates SNP lists per selection-enriched and unenriched pathway
as determined by GSEA}
\usage{
writePathFiles(
  genotype_file,
  results_file,
  gseaStat_file,
  snp2gene_file,
  fam_file,
  EM_group_file,
  ENRICH_NES = 3,
  UNENRICH_NES = 0.1,
  enrich_folder,
  enrichEM_folder,
  unenrich_folder
)
}
\arguments{
\item{genotype_file}{(char) path to file with SNP genotype data (PLINK format).}

\item{results_file}{(char) path to files with GSEA results.
Strutured to compare results of two population analyses
i.e., CEU vs. YRI and CEU vs. LWK.}

\item{gseaStat_file}{(char) path to GSEA statistics file.}

\item{snp2gene_file}{(char) path to SNP-gene mapping file.}

\item{fam_file}{(char) path to PLINK population coded fam file.}

\item{EM_group_file}{(char) file path to write pathway groupings as determined
by EnrichmentMap.}

\item{ENRICH_NES}{(integer) NES cutoff to select validated selection-enriched
pathways (default=3).}

\item{UNENRICH_NES}{(integer) NES cutoff to select unenriched pathways
(default=0.1).}

\item{enrich_folder}{(char) path to directory to store output files
(PLINK files per selection-enriched gene set).}

\item{enrichEM_folder}{(char) path to directory to store output files
(PLINK files per selection-enriched pathway, grouped via AutoAnnotate).}

\item{unenrich_folder}{(char) path to directory to store output files
(PLINK files per unenriched gene set).}
}
\value{
none
}
\description{
Generates SNP lists per selection-enriched and unenriched pathway
as determined by GSEA
}
