# [workdir]/ld/res...
# https://moderndata.plot.ly/power-curves-r-plotly-ggplot2/

library(pwr) # for power calcs
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation
library(ggplot2) # for plotting power curves
library(RColorBrewer)

hc <- readRDS("hc.diff.pairs.rds")
lc <- readRDS("lc.diff.pairs.rds")

hc_ixns <- as.numeric(table(hc$pathway))
lc_ixns <- nrow(lc)

powerTest <- function(sig, alt) {
  for (i in seq(0,1, length.out = 200)) {
    for (j in 1:length(unique(hc$pathway))) {
      pwrt[[j]] <- pwr.2p2n.test(h=i, n1=hc_ixns[j], n2=lc_ixns,
                                 sig.level=as.numeric(sprintf("%g", sig)),
                                 power=NULL,
                                 alternative=sprintf("%s", alt))

    }
    ptab <- rbind(ptab, cbind(pwrt[[1]]$h, pwrt[[1]]$power,
                              pwrt[[2]]$h, pwrt[[2]]$power,
                              pwrt[[3]]$h, pwrt[[3]]$power,
                              pwrt[[4]]$h, pwrt[[4]]$power,
                              pwrt[[5]]$h, pwrt[[5]]$power,
                              pwrt[[6]]$h, pwrt[[6]]$power,
                              pwrt[[7]]$h, pwrt[[7]]$power,
                              pwrt[[8]]$h, pwrt[[8]]$power,
                              pwrt[[9]]$h, pwrt[[9]]$power,
                              pwrt[[10]]$h, pwrt[[10]]$power,
                              pwrt[[11]]$h, pwrt[[11]]$power,
                              pwrt[[12]]$h, pwrt[[12]]$power,
                              pwrt[[13]]$h, pwrt[[13]]$power,
                              pwrt[[14]]$h, pwrt[[14]]$power,
                              pwrt[[15]]$h, pwrt[[15]]$power,
                              pwrt[[16]]$h, pwrt[[16]]$power,
                              pwrt[[17]]$h, pwrt[[17]]$power,
                              pwrt[[18]]$h, pwrt[[18]]$power,
                              pwrt[[19]]$h, pwrt[[19]]$power,
                              pwrt[[20]]$h, pwrt[[20]]$power,
                              pwrt[[21]]$h, pwrt[[21]]$power,
                              pwrt[[22]]$h, pwrt[[22]]$power,
                              pwrt[[23]]$h, pwrt[[23]]$power,
                              pwrt[[24]]$h, pwrt[[24]]$power,
                              pwrt[[25]]$h, pwrt[[25]]$power,
                              pwrt[[26]]$h, pwrt[[26]]$power,
                              pwrt[[27]]$h, pwrt[[27]]$power,
                              pwrt[[28]]$h, pwrt[[28]]$power,
                              pwrt[[29]]$h, pwrt[[29]]$power,
                              pwrt[[30]]$h, pwrt[[30]]$power,
                              pwrt[[31]]$h, pwrt[[31]]$power,
                              pwrt[[32]]$h, pwrt[[32]]$power,
                              pwrt[[33]]$h, pwrt[[33]]$power,
                              pwrt[[34]]$h, pwrt[[34]]$power,
                              pwrt[[35]]$h, pwrt[[35]]$power,
                              pwrt[[36]]$h, pwrt[[36]]$power,
                              pwrt[[37]]$h, pwrt[[37]]$power,
                              pwrt[[38]]$h, pwrt[[38]]$power,
                              pwrt[[39]]$h, pwrt[[39]]$power,
                              pwrt[[40]]$h, pwrt[[40]]$power,
                              pwrt[[41]]$h, pwrt[[41]]$power,
                              pwrt[[42]]$h, pwrt[[42]]$power,
                              pwrt[[43]]$h, pwrt[[43]]$power,
                              pwrt[[44]]$h, pwrt[[44]]$power,
                              pwrt[[45]]$h, pwrt[[45]]$power,
                              pwrt[[46]]$h, pwrt[[46]]$power,
                              pwrt[[47]]$h, pwrt[[47]]$power
            ))

  }

  ptab <- cbind(seq_len(nrow(ptab)), ptab)

  colnames(ptab) <- c("id",
                      "1_effect.size", "1_power",
                      "2_effect.size", "2_power",
                      "3_effect.size", "3_power",
                      "4_effect.size", "4_power",
                      "5_effect.size", "5_power",
                      "6_effect.size", "6_power",
                      "7_effect.size", "7_power",
                      "8_effect.size", "8_power",
                      "9_effect.size", "9_power",
                      "10_effect.size", "10_power",
                      "11_effect.size", "11_power",
                      "12_effect.size", "12_power",
                      "13_effect.size", "13_power",
                      "14_effect.size", "14_power",
                      "15_effect.size", "15_power",
                      "16_effect.size", "16_power",
                      "17_effect.size", "17_power",
                      "18_effect.size", "18_power",
                      "19_effect.size", "19_power",
                      "20_effect.size", "20_power",
                      "21_effect.size", "21_power",
                      "22_effect.size", "22_power",
                      "23_effect.size", "23_power",
                      "24_effect.size", "24_power",
                      "25_effect.size", "25_power",
                      "26_effect.size", "26_power",
                      "27_effect.size", "27_power",
                      "28_effect.size", "28_power",
                      "29_effect.size", "29_power",
                      "30_effect.size", "30_power",
                      "31_effect.size", "31_power",
                      "32_effect.size", "32_power",
                      "33_effect.size", "33_power",
                      "34_effect.size", "34_power",
                      "35_effect.size", "35_power",
                      "36_effect.size", "36_power",
                      "37_effect.size", "37_power",
                      "38_effect.size", "38_power",
                      "39_effect.size", "39_power",
                      "40_effect.size", "40_power",
                      "41_effect.size", "41_power",
                      "42_effect.size", "42_power",
                      "43_effect.size", "43_power",
                      "44_effect.size", "44_power",
                      "45_effect.size", "45_power",
                      "46_effect.size", "46_power",
                      "47_effect.size", "47_power")

    # get data into right format for ggplot2
    dat <<- ptab %>%
      as.data.frame() %>%
      gather(key = name, value = val, 2:ncol(ptab)) %>%
      separate(col = name, into = c("group", "var"), sep = "\\_") %>%
      spread(key = var, value = val)

    # factor group
    dat$group <- factor(dat$group, levels=1:47)
}

# Bonferroni-corrected sig level
pwrt <- list()
ptab <- cbind(NULL, NULL)
dat <- c()
alt <- "two.sided"
sig <- 0.05/length(unique(hc$pathway))

powerTest(sig=sig, alt=alt)
dat_bfcor <- dat
dat_bfcor$sig <- "Bonferroni"

# Nominal sig level
pwrt <- list()
ptab <- cbind(NULL, NULL)
dat <- c()
alt <- "two.sided"
sig <- 0.05

powerTest(sig=sig, alt=alt)
dat_nom <- dat
dat_nom$sig <- "Statistical"

both <- rbind(dat_nom, dat_bfcor)

# plot
cols <- colorRampPalette(brewer.pal(8, "Accent"))
npal <- cols(length(unique(hc$pathway)))

title <- paste("Power calculation comparing the proportion of\ninter-chromosomal",
               "associations per enriched\npathway against the nonenriched")
p <- ggplot(both, aes(x=effect.size, y=power, color=group)) +
     facet_grid(sig ~ .) +
     geom_line(size=1) +
     ggtitle(title) +
     labs(x="Effect size", y="Power") +
     theme(axis.text=element_text(size=14),
           axis.title=element_text(size=14),
           legend.position="none") +
    # geom_vline(xintercept = .54, linetype = 2) +
     geom_hline(yintercept = 0.80, linetype = 2) +
     scale_colour_manual(values=npal, guide=FALSE) +
     theme_Publication()

ggsave("power_2prop_test_sig-both_two-sided.png", p, width=7, height=8)
