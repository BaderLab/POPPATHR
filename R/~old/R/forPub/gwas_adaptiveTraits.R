# ~/GIN/data/geneProperties
gwas <- read.delim("nhgri_gwas_hits.txt", h=T, as.is=T)
toMatch <- paste(c("AIDS", "HIV", "hepatitis B", "hepatitis C", "tuberculosis",
                   "prion diseases", "malaria", "leprosy", "smallpox", "dengue",
                   "immune response", "anthrax", "Leishmaniasis", "influenza",
                   "human papillomavirus", "Helicobacter Pylori", "Meningococcal",
                   "cytomegalovirus", "skin pigment", "milk allergy", "height",
                   "hemoglobin level", "lipid metabolism", "altitude",
                   "hair morphology", "arsenic", "growth factor"),
            collapse="|")

test <- gwas[grep(toMatch, gwas$DISEASE.TRAIT, ignore.case=T),]
blah <- test[,c("PUBMEDID", "DATE", "STUDY", "DISEASE.TRAIT", "SNPS", "INITIAL.SAMPLE.SIZE")]
blah2 <- aggregate(.~PUBMEDID, blah, function(x) toString(unique(x)))
nrow(blah2)

snp_num <- lapply(strsplit(blah2$SNPS, ", "), length)
blah2$SNPS <- unlist(snp_num)

blah2 <- blah2[order(blah2$DISEASE.TRAIT),]
