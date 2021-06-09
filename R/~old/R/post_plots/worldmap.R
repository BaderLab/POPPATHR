# world map
# https://www.r-graph-gallery.com/how-to-draw-connecting-routes-on-map-with-r-and-great-circles/
library(tidyverse)
library(maps)
library(geosphere)

# draw an empty map
par(mar=c(0,0,0,0))
map('world', col="#C9CACC", fill=TRUE, bg="white", lwd=0.05,
    mar=rep(0,4), border=0, ylim=c(-60,85) )

# add city lat/long coordinates
# c([+east/-west], [+north/-south])
# HapMap populations
ASW = c(-110, 34) # southwest usa (african)
CEU = c(-111, 39) # utah (european)
CHB = c(116, 40) # beijing
CHD = c(-105, 40) # denver (chinese)
GIH = c(-95, 28) # houston (indian)
JPT = c(138, 36) # japan
LWK = c(35, 0) # webuye
MXL = c(-118, 34) # los angeles (mexican)
MKK = c(35, -1) # kinyawa (masaai)
TSI = c(11, 44) # toscani
YRI = c(4, 7) # ibadan (yoruba)

data <- rbind(ASW, CEU, CHB, CHD, GIH, JPT, LWK,
              MXL, MKK, TSI, YRI) %>% as.data.frame()
colnames(data) <- c("long", "lat")

# add points to map
points(x=data$long, y=data$lat, col="#662506", cex=1, pch=20)
#points(x=data$long, y=data$lat, col="#47C9EF", cex=1, pch=20)

# 1000 Genome populations
CHS = c(113, 22) # Southern Han Chinese
CDX = c(100, 22) # Chinese Dai in Xishuangbanna, China
KHV	= c(106, 10) # Kinh in Ho Chi Minh City, Vietnam
FIN	= c(25, 61) # Finnish in Finland
GBR	= c(-4, 56) # British in England and Scotland
IBS	= c(-3, 40) # Iberian Population in Spain
GWD	= c(-15, 13) # Gambian in Western Divisions in the Gambia
MSL	= c(-11, 8) # Mende in Sierra Leone
ESN	= c(8, 9) # Esan in Nigeria
ACB	= c(-59, 13) # African Caribbeans in Barbados
PUR	= c(-66, 18) # Puerto Ricans from Puerto Rico
CLM	= c(-75, 6) # Colombians from Medellin, Colombia
PEL	= c(-77, -12) # Peruvians from Lima, Peru
PJL	= c(74, 31) # Punjabi from Lahore, Pakistan
BEB	= c(90, 23) # Bengali from Bangladesh
STU	= c(-3, 55) # Sri Lankan Tamil from the UK
ITU	= c(-4, 55) # Indian Telugu from the UK

data2 <- rbind(CHS, CDX, KHV, FIN, GBR, IBS, GWD, MSL, ESN, ACB, PUR,
               CLM, PEL, PJL, BEB, STU, ITU) %>% as.data.frame()
colnames(data2) <- c("long", "lat")

points(x=data2$long, y=data2$lat, col="#7fc97f", cex=1, pch=20)

# drawing connections between population points
# get all pairs of coordinates
all_pairs <- cbind(t(combn(data$long, 2)),
             t(combn(data$lat, 2))) %>% as.data.frame()
colnames(all_pairs) <- c("long1", "long2", "lat1", "lat2")

# function that keeps the good part of the great circle, by Jeff Leek
getGreatCircle = function(userLL,relationLL){
  tmpCircle = greatCircle(userLL,relationLL, n=200)
  start = which.min(abs(tmpCircle[,1] - data.frame(userLL)[1,1]))
  end = which.min(abs(tmpCircle[,1] - relationLL[1]))
  greatC = tmpCircle[start:end,]
  return(greatC)
}

#plot connections
pdf("map.pdf")
par(mar=c(0,0,0,0))
map('world', col="#C9CACC", fill=TRUE, bg="white", lwd=0.05,
    mar=rep(0,4), border=0, ylim=c(-80,80) )

great <- getGreatCircle(ASW, CEU); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, CHB); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, CHD); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, GIH); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, JPT); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, LWK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, MXL); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(ASW, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, CHB); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, CHD); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, GIH); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, JPT); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, LWK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, MXL); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CEU, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, CHD); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, GIH); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, JPT); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, LWK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, MXL); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHB, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHD, GIH); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHD, JPT); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHD, LWK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHD, MXL); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHD, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHD, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(CHD, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(GIH, JPT); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(GIH, LWK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(GIH, MXL); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(GIH, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(GIH, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(GIH, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(JPT, LWK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(JPT, MXL); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(JPT, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(JPT, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(JPT, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(LWK, MXL); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(LWK, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(LWK, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(LWK, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(MXL, MKK); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(MXL, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(MXL, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(MKK, TSI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(MKK, YRI); lines(great, col="#47C9EF", lwd=1)
great <- getGreatCircle(TSI, YRI); lines(great, col="#47C9EF", lwd=1)

# draw points
points(x=data$long, y=data$lat, col="#FCB83F", cex=1, pch=20)
dev.off()
