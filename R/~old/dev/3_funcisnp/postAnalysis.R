# Post analysis script after running main FunciSNP programs
#pathGroup <- "high_conf"
dataDir <- paste0("/media/catherine/DATAPART1/Data/PopulationPathways",
                  "/methods/3_funciSNP/outfiles/GTEx_eQTLs",
                  "/out_2017-03-22_thyroid_eQTL_rand")
allFStableF <- sprintf("%s/all_tables/*.txt", dataDir)

#----------------------
suppressMessages({
    require(FunciSNP)
    require(dplyr)
    require(xlsx) })
#----------------------

# Combine all output tables to one table
# use awk to concatenate files while keeping single header
cat("\nCombining all FS tables...")
str1 <- "awk 'NR==1 {header=$_} FNR==1 && NR!=1"
str2 <- "{ $_ ~ $header getline; } {print}'"
str3 <- sprintf("%s > %s/FStable_all.txt", allFStableF, dataDir)
cmd <- sprintf("%s %s %s", str1, str2, str3)
system(cmd)
cat(sprintf(" file written to %s/FStable_all.txt\n", dataDir))

# Run post analysis
masterTableF <- sprintf("%s/FStable_all.txt", dataDir)
dat <- read.table(masterTableF, h=T, as.is=T) # use fread for larger tables
                                              # from data.table package
# Summary statistics
# Create FStable for high Rsq SNPs
rsq.filter = 0.5
FunciSNPtable(dat, rsq=rsq.filter)
highRsq <- filter(dat, R.squared >= rsq.filter)
write.table(highRsq,
            file=sprintf("%s/FStable_all_rsq_%g.txt", dataDir, rsq.filter),
            col=T, row=F, quote=F, sep="\t")

# Write summary table with tag SNP sums per biofeature
sumtable <- addmargins(table(highRsq$tag.snp.id, highRsq$bio.feature))
write.table(sumtable,
            file=sprintf("%s/FStable_summary_rsq_%g.txt", dataDir, rsq.filter),
            col=T, row=T, quote=F, sep="\t")

# Write to Excel format
df <- read.table(file=sprintf("%s/FStable_summary_rsq_%g.txt",
                 dataDir, rsq.filter), h=T, as.is=T)
write.xlsx(df,
           file=sprintf("%s/FStable_summary_rsq_%g.xlsx", dataDir, rsq.filter),
           col=T, row=T)
