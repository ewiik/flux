## testing schmidt stability measurement for buffalo pound for 2014
## the function requires the following units: water temperature in degrees C, depths (in m)
##    cross sectional areas (m^2), salinity (optional) in Practical Salinity Scale units
## Schmidt stability is in units J/m^2

## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  source("../scripts/diel-buffalo.R")}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- bdat[order(bdat$datetime),]

bdat5 <- readRDS('../data/private/bpbuoy2015-mod.rds')
bdat5 <- bdat5[order(bdat5$datetime),]

if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  stop("get hypsometry from Emma")}
hyp <- read.csv('../data/private/BuffaloPound_Hypso-1.csv')

## load packages
stopifnot(require("rLakeAnalyzer"))
stopifnot(require("reshape2"))
stopifnot(require("ggplot2"))
stopifnot(require("viridis"))
stopifnot(require("gridExtra"))

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=18, base_family = 'Arial') +
  theme(legend.position='top')

## melt vectors for the temperature profile : 2014
btemp <- melt(bdat, id.vars = "datetime", measure.vars = c("temp1", "temp2", "temp4",
                                                           "temp5", "temp6"))
btemp$depth <- rep(0) #45, 77, 123, 218 + reg 82 depth
btemp$depth[grep("1", btemp$variable)] <- 0.82
btemp$depth[grep("2", btemp$variable)] <- 2.94
btemp$depth[grep("4", btemp$variable)] <- 0.77
btemp$depth[grep("5", btemp$variable)] <- 1.23
btemp$depth[grep("6", btemp$variable)] <- 2.18

## melt vectors for the temperature profile : 2015 (>3000 NAs for temp1, not taking it)
btemp5 <- melt(bdat5, id.vars = "datetime", measure.vars = c("temp2", "temp3", "temp4",
                                                           "temp5", "temp6", "temp7"))
btemp5$depth <- rep(0) #0.45, 0.77, 1.23, 2.18, 3.18 + reg 82, 2.92 depth
btemp5$depth[grep("2", btemp5$variable)] <- 2.94
btemp5$depth[grep("3", btemp5$variable)] <- 0.45
btemp5$depth[grep("4", btemp5$variable)] <- 0.77
btemp5$depth[grep("5", btemp5$variable)] <- 1.23
btemp5$depth[grep("6", btemp5$variable)] <- 2.18
btemp5$depth[grep("7", btemp5$variable)] <- 3.18

## schmidt function does not calculate anything if any one value is NA. Let's remove
##    these from btemp
btemp <- na.omit(btemp)
btemp5 <- na.omit(btemp5)

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

## split to lists with vectors of each variable: 2014
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

## split to lists with vectors of each variable: 2015
testing5 <- with(btemp5, split(btemp5, list(datetime)))
templist5 <- lapply(testing, '[[', 3)
tempdepthlist5 <- lapply(testing, '[[', 4)
sallist5 <- lapply(tempdepthlist, function(x) x*0)

bthA <- hyp$revcum
bthD <- hyp$depth
arealist5 <- vector(mode="list", length=nrow(bdat5))
arealist5[1:length(arealist)] <- list(bthA)
depthlist5 <- vector(mode="list", length=nrow(bdat5))
depthlist5[1:length(depthlist)] <- list(bthD)

## apply schmidt function to each date point:2014
stablist <- lapply(mapply(schmidt.stability, wtr=templist, depths=tempdepthlist,
                                  bthA=arealist, bthD=depthlist, sal=sallist), '[[', 1)
stablist <- as.data.frame(do.call(rbind, stablist))
rownames(stablist) <- NULL
names(stablist) <- "stability"
stablist$datetime <- bdat$datetime

## apply schmidt function to each date point:2015
stablist5 <- lapply(mapply(schmidt.stability, wtr=templist5, depths=tempdepthlist5,
                          bthA=arealist, bthD=depthlist5, sal=sallist5), '[[', 1)
stablist5 <- as.data.frame(do.call(rbind, stablist5))
rownames(stablist5) <- NULL
names(stablist5) <- "stability"
stablist5$datetime <- bdat5$datetime

## save for modeling
saveRDS(stablist, "../data/private/bp-stability.rds")
saveRDS(stablist5, "../data/private/bp-stability2015.rds")

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

