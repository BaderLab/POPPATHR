library(plyr)
library(ggplot2)

dataDir <- "/media/catherine/DATAPART1/Data/PopulationPathways/methods/2_snpstats/data_rds/pnc"
ceuDir <- sprintf("%s/ceu", dataDir)
aswDir <- sprintf("%s/asw", dataDir)

#load in data
hc.diff.r2.mean.pop1 <- readRDS(sprintf("%s/hc_diff_r2_ceu.rds", ceuDir))
hc.diff.dp.mean.pop1 <- readRDS(sprintf("%s/hc_diff_dprime_ceu.rds", ceuDir))
hc.diff.r2.mean.pop2 <- readRDS(sprintf("%s/hc_diff_r2_asw.rds", aswDir))
hc.diff.dp.mean.pop2 <- readRDS(sprintf("%s/hc_diff_dprime_asw.rds", aswDir))

null.diff.r2.mean.pop1 <- readRDS(sprintf("%s/null_diff_r2_500x400_ceu.rds", ceuDir))
null.diff.dp.mean.pop1 <- readRDS(sprintf("%s/null_diff_dprime_500x400_ceu.rds", ceuDir))
null.diff.r2.mean.pop2 <- readRDS(sprintf("%s/null_diff_r2_500x400_asw.rds", aswDir))
null.diff.dp.mean.pop2 <- readRDS(sprintf("%s/null_diff_dprime_500x400_asw.rds", aswDir))

real.pop1 <- data.frame(R.squared=hc.diff.r2.mean.pop1,
                        D.prime=hc.diff.dp.mean.pop1,
                        pathway.group="highconf")

null.pop1 <- data.frame(R.squared=null.diff.r2.mean.pop1,
                        D.prime=null.diff.dp.mean.pop1,
                        pathway.group="rand")

dist.dat.pop1 <- rbind(real.pop1, null.pop1)
mean.dat.pop1 <- ddply(dist.dat.pop1, "pathway.group", summarise,
                       R.squared.mean=mean(R.squared),
                       D.prime.mean=mean(D.prime))
####
real.pop2 <- data.frame(R.squared=hc.diff.r2.mean.pop2,
                       D.prime=hc.diff.dp.mean.pop2,
                       pathway.group="highconf")

null.pop2 <- data.frame(R.squared=null.diff.r2.mean.pop2,
                       D.prime=null.diff.dp.mean.pop2,
                       pathway.group="rand")

dist.dat.pop2 <- rbind(real.pop2, null.pop1)
mean.dat.pop2 <- ddply(dist.dat.pop2, "pathway.group", summarise,
                      R.squared.mean=mean(R.squared),
                      D.prime.mean=mean(D.prime))

##################################
title <- paste("Degree of co-selection per interchromosomal SNP-SNP",
               "\ninteraction within the high-confidence pathways",
               "\nvs. randomly selected SNP groups")
r2.axis.title <- bquote("Pairwise LD value (mean " *r^2*" per pathway)")
dp.axis.title <- bquote("Pairwise LD value (mean D' per pathway)")

## for use with stat_summary(fun.data=box.style); allows white median line to
## appear after colouring and filling boxplots
box.style <- function(x){
    return(c(y=median(x), ymin=median(x), ymax=median(x)))
}

## for use with stat.summary(fun.data=give.n); displays sample size (N)
## courtesy of Bangyou at Stack Overflow
give.n <- function(x){
    return(c(y=median(x)*1.50, label=length(x)))
    # experiment with the multiplier to find the perfect position
}

## calculate empirical p by quantifying all permuted mean r2 values greater
## than the real mean r2 values, divided by the total number of replicates
calc.p <- function(null, real){
  return( (length(which(null > mean(real)))+1) / (length(null)+1) )
}

###################################
ggplot(dist.dat.pop1, aes(x=R.squared, colour=pathway.group,
                          fill=pathway.group)) +
    geom_density(alpha=0.3) +
    #geom_histogram(bins=20, alpha=0.5, position="identity") +
    geom_vline(data=mean.dat.pop1, aes(xintercept=R.squared.mean,
                                       colour=pathway.group),
               linetype="dashed", size=1) +
    scale_x_continuous(r2.axis.title) +
    scale_y_continuous("Density") +
    ggtitle(title) +
    theme_set(theme_minimal()) +
    theme(plot.title=element_text(hjust=0.5),
          text=element_text(size=17),
          legend.position="top",
          legend.title=element_blank(),
          panel.grid.major.x=element_blank())
ggsave("dist_nullvreal_r2_ceu.png", width=8, height=7)

perm.p = calc.p(null.diff.r2.mean.pop1, hc.diff.r2.mean.pop1)
cat(sprintf("P-value for permuted sample vs. real test statistic= %g\n",
            perm.p))

###################################
ggplot(dist.dat.pop1, aes(x=pathway.group, y=R.squared)) +
      geom_boxplot(outlier.colour=NULL,
                   aes(colour=pathway.group, fill=pathway.group)) +
      stat_summary(geom="crossbar", width=0.65, fatten=0, color="white",
                   fun.data=box.style) +
      scale_y_continuous(r2.axis.title) +
      ggtitle(title) +
      theme_set(theme_minimal()) +
      theme(plot.title=element_text(hjust=0.5),
            text=element_text(size=17),
            legend.position="top",
            legend.title=element_blank(),
            panel.grid.major.x=element_blank(),
            axis.title.x=element_blank()) +
      scale_x_discrete(labels=paste("N=", table(dist.dat.pop1$pathway.group),
                                    sep=""))
  ggsave("boxplot_nullvsreal_r2_ceu.png", width=8, height=7.5)
