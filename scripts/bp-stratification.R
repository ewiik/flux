## testing schmidt stability measurement for buffalo pound for 2014
## the function requires the following units: water temperature in degrees C, depths (in m) 
##    cross sectional areas (m^2), salinity (optional) in Practical Salinity Scale units
## Schmidt stability is in units J/m^2

## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  source("../scripts/diel-buffalo.R")}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- bdat[order(bdat$datetime),]

if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  stop("get hypsometry from Emma")}
hyp <- read.csv('../data/private/BuffaloPound_Hypso-1.csv')

## load packages
library("rLakeAnalyzer")
library("reshape2")
library("ggplot2")
library("viridis")
library("gridExtra")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=18, base_family = 'Arial') +
  theme(legend.position='top')

## melt vectors for the temperature profile
btemp <- melt(bdat, id.vars = "datetime", measure.vars = c("temp1", "temp2", "temp4",
                                                           "temp5", "temp6"))
btemp$depth <- rep(0) #45, 77, 123, 218 + reg 82 depth
btemp$depth[grep("1", btemp$variable)] <- 0.82
btemp$depth[grep("2", btemp$variable)] <- 2.94
btemp$depth[grep("4", btemp$variable)] <- 0.77
btemp$depth[grep("5", btemp$variable)] <- 1.23
btemp$depth[grep("6", btemp$variable)] <- 2.18

btemp$value[which(btemp$value < 0)] <- NA
btemp$value[which(btemp$value > 23 & btemp$variable == "temp6")] <- NA 
# no idea why this is happening in temp6!

## schmidt function does not calculate anything if any one value is NA. Let's remove 
##    these from btemp
btemp <- na.omit(btemp)

## arange hyp into format required: in examples in pdf,
##    bthA    <- c(10000,8900,5000,3500,2000,1000,300,10)
##    bthD    <- c(0,1,2,3,4,5,6,7
##    which implies that i need to start with 0 has full area, etc., etc.
revcum <- function(area) {
  area <- area
  target <- vector('numeric', length=length(area))
  for (i in 1:length(area)) {
    target[i] <- sum(area[i:length(area)])
  }
  target
}
hyp <- hyp[,c('Bins_Depth.m.',"BinFrequency")]
names(hyp) <- c("depth", "area")
hyp <- na.omit(hyp)
hyp <- hyp[-which(hyp$area == 0),]

hyp$revcum <- with(hyp, revcum(area))

## split to lists with vectors of each variable
testing <- with(btemp, split(btemp, list(datetime))) 
templist <- lapply(testing, '[[', 3)
tempdepthlist <- lapply(testing, '[[', 4)
sallist <- lapply(tempdepthlist, function(x) x*0)

bthA <- hyp$revcum
bthD <- hyp$depth
arealist <- vector(mode="list", length=nrow(bdat))
arealist[1:length(arealist)] <- list(bthA)
depthlist <- vector(mode="list", length=nrow(bdat))
depthlist[1:length(depthlist)] <- list(bthD)

## apply schmidt function to each date point
stablist <- lapply(mapply(schmidt.stability, wtr=templist, depths=tempdepthlist, 
                                  bthA=arealist, bthD=depthlist, sal=sallist), '[[', 1)
stablist <- as.data.frame(do.call(rbind, stablist))
rownames(stablist) <- NULL
names(stablist) <- "stability"
stablist$datetime <- bdat$datetime

## merge with bdat for plotting
tester <- merge(stablist, bdat)

## compare with temperature profile
highstab <- which(tester$stability >= 5)
highdates <- tester$datetime[highstab]
highdates <- data.frame(x=highdates, y=rep(16))
stabplot <- ggplot(tester, aes(y=stability, x=datetime)) +
  papertheme +
  geom_line()
tempplot <- ggplot(btemp, aes(y=value, x=datetime, col=depth)) +
  papertheme +
  geom_point(size=0.3) +
  geom_point(inherit.aes = FALSE, data=highdates, aes(x=x, y=y), size=0.3, show.legend = FALSE) +
  scale_color_viridis(discrete = FALSE, begin= 0.1, end=0.9, alpha=0.4,
                      option = 'viridis', direction = -1, name="Depth (m)")  +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank()) +
  ylab("Temperature")
grid.arrange(tempplot, stabplot)  

