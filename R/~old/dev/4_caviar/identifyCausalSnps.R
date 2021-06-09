#' Quantify the probability of causal variants using CAVIAR
#' (http://genetics.cs.ucla.edu/caviar/)


# LD matrix
./plink --bfile pnc_subset_tagSnps --r2 square spaces --ld-window-r2 0.2
        --allow-no-sex --out tagSnps_LDmatrix

#square - matrix format (needed for caviar)
#space - space delimited (as per example format)

#get MAF values for tag SNPs
plink <- read.table("...bim", h=F, as.is=T)
maf <- read.table("markerMAF", h=T, as.is=T)
colnames(plink)[2] <- "Marker"
plink <- as.data.frame(plink[2])
topsnps <- merge(plink, maf, by="Marker")
write.table(topsnps, "topSnps_MAFscores.txt", col=F, row=F, quote=F, sep="\t")

#run caviar
./CAVIAR -z tagSnps_MAF.txt -l tagSnps_LDmatrix.ld -o test
