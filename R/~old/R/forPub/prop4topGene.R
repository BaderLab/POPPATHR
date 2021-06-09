
c1 <- paste(c('CXCL2', 'KRAS', 'VTI1A', 'PTPN6', 'CD86', 'PEX14'), collapse="|")
c2 <- paste(c('FGF1', 'GRB2', 'ACVR1C', 'STAT3', 'PLCB1'), collapse="|")
c3 <- paste(c('MTOR', 'BRCA1', 'RNF145', 'LDLR'), collapse="|")
c4 <- paste(c('IFI16', 'RHO', 'NOC2L', 'RAD9A', 'TP53'), collapse="|")
c5 <- paste(c('NOTCH2', 'FOXP1', 'IRF4', 'LAMC1'), collapse="|")
c6 <- paste(c('LCOR', 'ESR1', 'PHB2', 'COL18A1'), collapse="|")
c7 <- paste(c('BAHD1', 'HDAC5', 'SIRT1', 'HIST1H4H'), collapse="|")

gene_list <- list(c1, c2, c3, c4, c5, c6, c7)

countDrugs <- function(genes) {
  blah <- as.character(all_prop[grep(genes, all_prop$gene),]$dgidb_drug_name)
  blah <- strsplit(blah, ", ")
  blah <- unique(unlist(blah))
  blah <- na.omit(blah)
  length(blah)
}

mapply(countDrugs, gene_list)
