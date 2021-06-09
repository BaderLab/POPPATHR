# [workdir]/ld/res...
# https://moderndata.plot.ly/power-curves-r-plotly-ggplot2/

library(pwr) # for power calcs
library(dplyr) # for data manipulation
library(tidyr) # for data manipulation
library(ggplot2) # for plotting power curves
library(RColorBrewer)

hc <- readRDS("hc.diff.pairs.rds")

# Generate power calculations
powerTest <- function(sig, alt) {
  for (i in seq(0,1, length.out = 200)){
    pwrt1 <- pwr.2p2n.test(h = i, n1 = 4766, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt2 <- pwr.2p2n.test(h = i, n1 = 10388, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt3 <- pwr.2p2n.test(h = i, n1 = 11368, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt4 <- pwr.2p2n.test(h = i, n1 = 932, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt5 <- pwr.2p2n.test(h = i, n1 = 2644, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt6 <- pwr.2p2n.test(h = i, n1 = 4700, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt7 <- pwr.2p2n.test(h = i, n1 = 15168, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt8 <- pwr.2p2n.test(h = i, n1 = 4541, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt9 <- pwr.2p2n.test(h = i, n1 = 2549, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt10 <- pwr.2p2n.test(h = i, n1 = 2420, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt11 <- pwr.2p2n.test(h = i, n1 = 8707, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt12 <- pwr.2p2n.test(h = i, n1 = 9581, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt13 <- pwr.2p2n.test(h = i, n1 = 2438, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt14 <- pwr.2p2n.test(h = i, n1 = 3811, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt15 <- pwr.2p2n.test(h = i, n1 = 15087, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt16 <- pwr.2p2n.test(h = i, n1 = 4433, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt17 <- pwr.2p2n.test(h = i, n1 = 12579, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt18 <- pwr.2p2n.test(h = i, n1 = 13285, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt19 <- pwr.2p2n.test(h = i, n1 = 12602, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))
    pwrt20 <- pwr.2p2n.test(h = i, n1 = 50, n2 = 18437,
                          sig.level = as.numeric(sprintf("%g", sig)),
                          power = NULL,
                          alternative=sprintf("%s", alt))

    ptab <<- rbind(ptab, cbind(pwrt1$h, pwrt1$power,
                              pwrt2$h, pwrt2$power,
                              pwrt3$h, pwrt3$power,
                              pwrt4$h, pwrt4$power,
                              pwrt5$h, pwrt5$power,
                              pwrt6$h, pwrt6$power,
                              pwrt7$h, pwrt7$power,
                              pwrt8$h, pwrt8$power,
                              pwrt9$h, pwrt9$power,
                              pwrt10$h, pwrt10$power,
                              pwrt11$h, pwrt11$power,
                              pwrt12$h, pwrt12$power,
                              pwrt13$h, pwrt13$power,
                              pwrt14$h, pwrt14$power,
                              pwrt15$h, pwrt15$power,
                              pwrt16$h, pwrt16$power,
                              pwrt17$h, pwrt17$power,
                              pwrt18$h, pwrt18$power,
                              pwrt19$h, pwrt19$power,
                              pwrt20$h, pwrt20$power
                    ))
          }
  ptab <<- cbind(seq_len(nrow(ptab)), ptab)

  colnames(ptab) <- c("id",
                      "n1=4766, n2=18437.effect size","n1=4766, n2=18437.power",
                      "n1=10388, n2=18437.effect size","n1=10388, n2=18437.power",
                      "n1=11368, n2=18437.effect size","n1=11368, n2=18437.power",
                      "n1=932, n2=18437.effect size","n1=932, n2=18437.power",
                      "n1=2644, n2=18437.effect size","n1=2644, n2=18437.power",
                      "n1=4700, n2=18437.effect size","n1=4700, n2=18437.power",
                      "n1=15168, n2=18437.effect size","n1=15168, n2=18437.power",
                      "n1=4541, n2=18437.effect size","n1=4541, n2=18437.power",
                      "n1=2549, n2=18437.effect size","n1=2549, n2=18437.power",
                      "n1=2420, n2=18437.effect size","n1=2420, n2=18437.power",
                      "n1=8707, n2=18437.effect size","n1=8707, n2=18437.power",
                      "n1=9581, n2=18437.effect size","n1=9581, n2=18437.power",
                      "n1=2438, n2=18437.effect size","n1=2438, n2=18437.power",
                      "n1=3811, n2=18437.effect size","n1=3811, n2=18437.power",
                      "n1=15087, n2=18437.effect size","n1=15087, n2=18437.power",
                      "n1=4433, n2=18437.effect size","n1=4433, n2=18437.power",
                      "n1=12579, n2=18437.effect size","n1=12579, n2=18437.power",
                      "n1=13285, n2=18437.effect size","n1=13285, n2=18437.power",
                      "n1=12602, n2=18437.effect size","n1=12602, n2=18437.power",
                      "n1=50, n2=50.effect size", "n1=50, n2=50.power")

  # get data into right format for ggplot2
  dat <<- ptab %>%
    as.data.frame() %>%
    gather(key = name, value = val, 2:ncol(ptab)) %>%
    separate(col = name, into = c("group", "var"), sep = "\\.") %>%
    spread(key = var, value = val)

  # factor group
  dat$group <- factor(dat$group,
                  levels = c("n1=4766, n2=18437",
                             "n1=10388, n2=18437",
                             "n1=11368, n2=18437",
                             "n1=932, n2=18437",
                             "n1=2644, n2=18437",
                             "n1=4700, n2=18437",
                             "n1=15168, n2=18437",
                             "n1=4541, n2=18437",
                             "n1=2549, n2=18437",
                             "n1=2420, n2=18437",
                             "n1=8707, n2=18437",
                             "n1=9581, n2=18437",
                             "n1=2438, n2=18437",
                             "n1=3811, n2=18437",
                             "n1=15087, n2=18437",
                             "n1=4433, n2=18437",
                             "n1=12579, n2=18437",
                             "n1=13285, n2=18437",
                             "n1=12602, n2=18437",
                             "n1=50, n2=18437"))

}

# Bonferroni-corrected sig level
ptab <- cbind(NULL, NULL)
dat <- c()
alt <- "two.sided"
sig <- 0.05/length(unique(hc$pathway))

powerTest(sig=sig, alt=alt)
dat_bfcor <- dat
dat_bfcor$sig <- "Bonferroni"

# Nominal sig level
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
npal <- cols(length(unique(hc$pathway))+1)

title <- paste("Power calculation comparing the proportion of inter-chromosomal\n",
               "interactions per ancestry-enriched pathway against the set\n",
               "of nonenriched pathways")
p <- ggplot(both, aes(x = `effect size`, y = power, color = group)) +
     facet_grid(sig ~ .) +
     geom_line(size=2) +
     ggtitle(title) +
     labs(x="Effect size", y="Power") +
     theme(axis.text=element_text(size=14),
           axis.title=element_text(size=14),
           legend.text=element_text(size=14)) +
    # geom_vline(xintercept = .54, linetype = 2) +
     geom_hline(yintercept = 0.80, linetype = 2) +
     scale_colour_manual(values=npal) +
     theme_Publication()

ggsave("power_test.png", p, width=8, height=8)
ggsave("power_2prop_test_sig-both_two-sided.png", p, width=8, height=8)
