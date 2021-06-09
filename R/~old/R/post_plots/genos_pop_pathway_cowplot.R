# [workdir]/ld
# get pop-based snp genotypes from snpMatrix data
require(snpStats)
require(reshape2)
require(ggplot2)
require(cowplot)

getGenos <- function(set) {
  bed <- list.files(path=sprintf("%s_snps", set), pattern="*.bed", full.names=T)
  bim <- list.files(path=sprintf("%s_snps", set), pattern="*.bim", full.names=T)
  fam <- list.files(path=sprintf("%s_snps", set), pattern="*.fam", full.names=T)

  freq <- c()
  for (i in 1:length(bed)) {
    dat <- read.plink(bed[i], bim[i], fam[i])

    # get population assignments
    pops <- dat$fam
    ceu <- pops$member[which(pops$affected == 1)]
    asw <- pops$member[which(pops$affected == 2)]

    # subset genotypes by population
    geno_mat <- as(dat$genotypes, "numeric")
    geno_mat_ceu <- geno_mat[ceu, ]
    geno_mat_asw <- geno_mat[asw, ]

    geno_ceu <- melt(geno_mat_ceu)
    geno_ceu <- na.omit(geno_ceu)
    geno_asw <- melt(geno_mat_asw)
    geno_asw <- na.omit(geno_asw)

    # re-format for plotting
    ceu_freq <- as.data.frame(table(geno_ceu$value)/nrow(geno_ceu))
    ceu_freq$Population <- "CEU"
    ceu_freq$pathway <- sprintf("pathway_%i", i)
    asw_freq <- as.data.frame(table(geno_asw$value)/nrow(geno_asw))
    asw_freq$Population <- "YRI"
    asw_freq$pathway <- sprintf("pathway_%i", i)

    freq[[i]] <- rbind(ceu_freq, asw_freq)
  }

  freq <<- do.call("rbind", freq)
    #freq_melt <- melt(freq, id.vars=c("Var1", "pop", "pathway"), measure.vars="Freq")
}

getGenos(set="hc")
freq_hc <- freq
freq_hc$set <- "Selection-enriched"

getGenos(set="lc")
freq_lc <- freq
freq_lc$set <- "Unenriched"

both <- rbind(freq_hc, freq_lc)

title <- "Distribution of pathway-level SNP genotype frequencies"
p <- ggplot(both, aes(x=as.factor(Var1), y=Freq, fill=Population)) +
      facet_wrap(~set, ncol=2) +
      geom_bar(stat="identity", position=position_dodge()) +
      ggtitle(title) +
      labs(x="Genotype", y="Relative percentage") +
      scale_fill_manual(values=c("#7D54A5", "#A6CEE3"))

save_plot("snp_genotypes_path_pop.tiff", p, base_width=8, base_height=4)
