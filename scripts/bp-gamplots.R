## plots of the bp diel data

## load packages
library('ggplot2')
library('mgcv')
library('cowplot')
library("extrafont")
library("viridis")
library("grid")
library("gridExtra")
library("reshape2")

source("../functions/geom_rug3.R")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

## read in models for both years
co2mod <- readRDS('../data/private/bpco2phmod.rds')

phmod <- readRDS('../data/private/BPphmod.rds')
phmod5 <- readRDS('../data/private/BPphmod2015.rds')

bptimemod <- readRDS("../data/private/bp-diel-timemod2014.rds")

bpco2mod <- readRDS("../data/private/bpmod.rds")
co2minus <- readRDS("../data/private/bpmodsimp.rds")

## read in data
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat5 <- readRDS('../data/private/bpbuoy2015-mod.rds')

if (!file.exists("../data/private/bp-stability2015.rds")) {
  source("./bp-stratification.R")
}
schmidt <- readRDS("../data/private/bp-stability.rds") # note this for 2014
bdat <- merge(bdat, schmidt)

schmidt5 <- readRDS("../data/private/bp-stability2015.rds") # note this for 2014
bdat5 <- merge(bdat5, schmidt5)

## create index of potential convective mixing, based on diff in T btw air and water
bdat$conv <- bdat$airtemp - bdat$temp4 # this is the topmost, 77cm one
bdat5$conv <- bdat5$airtemp - bdat5$temp4 # this is the topmost, 77cm one


#bdat$TimeofDay <- as.factor(bdat$TimeofDay)

## predict for all vars in co2minus
## ============================================================================================================
## limit prediction data frame to observed intervals
regsplit <- with(bdat, split(bdat, list(TimeofDay)))
regsplit <- regsplit[sapply(regsplit, function(x) dim(x)[1]) > 0] #remove empties for lakes R, E etc.
minmax <- function(df, colnames) {
  allmin <- as.data.frame(do.call(cbind, lapply(df[,colnames], min, na.rm = TRUE)))
  names(allmin) <- sapply(names(allmin), function(x) paste("min",x, sep = ""))
  allmax <- as.data.frame(do.call(cbind, lapply(df[,colnames], max, na.rm = TRUE)))
  names(allmax) <- sapply(names(allmax), function(x) paste("max",x, sep = ""))
  summ <- as.data.frame(cbind(allmin, allmax))
  summ <- cbind(summ, data.frame(TimeofDay = df$TimeofDay[1]))
  summ
}
minmaxes <- do.call(rbind, lapply(regsplit, minmax, colnames = c("chlrfu","bga1rfu", "airtemp", "stability", 
                                                                 "dailyrain", "conv", "turb","cdom",
                                                                 "windsp", "cond1")))
rownames(minmaxes) <- NULL


## ==================================================================================================
## 1. chl 
## ===================================================================================================
N <- 200
varWant <- c("bga1rfu", "airtemp", "stability", "dailyrain", "conv", "turb","cdom","windsp",
             "cond1")
lakeXbar <- with(bdat, do.call("rbind",
                                   lapply(split(bdat[, varWant], droplevels(TimeofDay)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, TimeofDay = factor(rownames(lakeXbar)))

co2.pdat <- with(droplevels(bdat),
                 data.frame(`chlrfu` = rep(seq(min(`chlrfu`, na.rm = TRUE),
                                                   max(`chlrfu`, na.rm = TRUE),
                                                   length = N),
                                               nlevels(TimeofDay)),
                            TimeofDay = rep(levels(TimeofDay), each = N)
                            ))
co2.pdat <- merge(co2.pdat, lakeXbar)

## let's see what the output is for 
co2.pred <- predict(co2minus, newdata = co2.pdat, type = "terms", se.fit = TRUE)
co2.test <- predict(co2minus, newdata = co2.pdat, type = "response", se.fit = TRUE)
co2.test2 <- predict(co2minus, newdata = co2.pdat, type = "link", se.fit = TRUE)
## ok this works: exp is the link function's return path. (see also mod$family$linkinv)
##    the following all produce the same
exp(head(co2.test2$fit))
head(co2.test$fit)
#exp(head(rowSums(co2.pred$fit)) + shiftco2) 
  
whichCols <- grep("chl", colnames(co2.pred$fit))
whichColsSE <- grep("chl", colnames(co2.pred$se.fit))

co2.resp <- cbind(co2.pdat, Fitted = rowSums(co2.pred$fit[, whichCols]), 
                  se.Fitted = rowSums(co2.pred$se.fit[, whichColsSE]))

co2.resp <- with(co2.resp, transform(co2.resp, Fittedplus = Fitted + se.Fitted))
co2.resp <- with(co2.resp, transform(co2.resp, Fittedminus = Fitted - se.Fitted))

## make into original limits
shiftco2 <- attr(predict(co2minus, newdata = co2.pdat, type = "iterms"), "constant")
co2.respnorm <- co2.resp
co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fitted = Fitted + shiftco2))
co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fittedplus = Fittedplus + shiftco2))
co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fittedminus = Fittedminus + shiftco2))
whichexp <- grep("Fitted", names(co2.respnorm))
co2.respnorm[,whichexp] <- exp(co2.respnorm[,whichexp])

co2.respnorm <- merge(co2.respnorm, minmaxes)
overs <- with(co2.respnorm, which(chlrfu < minchlrfu | chlrfu > maxchlrfu))
co2.respnorm <- co2.respnorm[-overs,]


chlplot <- ggplot(co2.respnorm, aes(x = chlrfu, y = Fitted, 
                                    colour = ifelse(TimeofDay == "Day", "Day", "Night"),
                                    lty=ifelse(TimeofDay == "Day", "Day", "Night"))) +
  papertheme + 
  #annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
             alpha = 0.25, fill='white') +  
  scale_linetype_manual(name='Time', values = c("solid", "dotdash")) +
  scale_colour_manual(name="Time", values = c("#b2abd2", "#5e3c99"))+ #"#b2abd2", "#5e3c99", "#e66101",
    #"#5e3c99", "#5e3c99"
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE, inherit.aes = FALSE) +
  #geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  #geom_vline(xintercept = meanpH, linetype = 'dotted') +
  geom_rug3(aes(x=chlrfu, y=co2corr), data = bdat, stat = "identity", position = "identity", 
            sides = "bl", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  theme(legend.position='top') +
  guides(colour=guide_legend(ncol=3, bycol =TRUE,title.position = 'left')) +
  xlab('Chlorophyll') + ylab(expression(paste(italic('p')*'CO'[2]~"(ppm)")))

## 2. bga1rfu
## ==================================================================================================
N <- 200
varWant <- c("chlrfu", "airtemp", "stability", "ODOrel1", "dailyrain", "conv", "turb","cdom","windsp",
             "cond1")
lakeXbar <- with(bdat, do.call("rbind",
                               lapply(split(bdat[, varWant], droplevels(TimeofDay)), 
                                      colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, TimeofDay = factor(rownames(lakeXbar)))

co2.pdat <- with(droplevels(bdat),
                 data.frame(`bga1rfu` = rep(seq(min(`bga1rfu`, na.rm = TRUE),
                                               20,
                                               length = N),
                                           nlevels(TimeofDay)),
                            TimeofDay = rep(levels(TimeofDay), each = N)
                 ))
co2.pdat <- merge(co2.pdat, lakeXbar)

## let's see what the output is for 
co2.pred <- predict(co2minus, newdata = co2.pdat, type = "terms", se.fit = TRUE)
co2.test <- predict(co2minus, newdata = co2.pdat, type = "response", se.fit = TRUE)
co2.test2 <- predict(co2minus, newdata = co2.pdat, type = "link", se.fit = TRUE)
## ok this works: exp is the link function's return path. (see also mod$family$linkinv)
##    the following all produce the same
exp(head(co2.test2$fit))
head(co2.test$fit)
exp(head(rowSums(co2.pred$fit)) + shiftco2) 

whichCols <- grep("bga", colnames(co2.pred$fit))
whichColsSE <- grep("bga", colnames(co2.pred$se.fit))

co2.resp <- cbind(co2.pdat, Fitted = rowSums(co2.pred$fit[, whichCols]), 
                  se.Fitted = rowSums(co2.pred$se.fit[, whichColsSE]))

co2.resp <- with(co2.resp, transform(co2.resp, Fittedplus = Fitted + se.Fitted))
co2.resp <- with(co2.resp, transform(co2.resp, Fittedminus = Fitted - se.Fitted))

## make into original limits
shiftco2 <- attr(predict(co2minus, newdata = co2.pdat, type = "iterms"), "constant")
co2.respnorm <- co2.resp
co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fitted = Fitted + shiftco2))
co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fittedplus = Fittedplus + shiftco2))
co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fittedminus = Fittedminus + shiftco2))
whichexp <- grep("Fitted", names(co2.respnorm))
co2.respnorm[,whichexp] <- exp(co2.respnorm[,whichexp])

cyanoplot <- ggplot(co2.respnorm, aes(x = bga1rfu, y = Fitted, 
                                    colour = ifelse(TimeofDay == "Day", "Day", "Night"),
                                    lty=ifelse(TimeofDay == "Day", "Day", "Night"))) +
  papertheme + 
  #annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25, fill='white') +  
  scale_linetype_manual(name='Time', values = c("solid", "dotdash")) +
  scale_colour_manual(name="Time", values = c("#b2abd2", "#5e3c99"))+ #"#b2abd2", "#5e3c99", "#e66101",
  #"#5e3c99", "#5e3c99"
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE, inherit.aes = FALSE) +
  #geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  #geom_vline(xintercept = meanpH, linetype = 'dotted') +
  geom_rug3(aes(x=chlrfu, y=co2corr), data = bdat, stat = "identity", position = "identity", 
            sides = "bl", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  theme(legend.position='top') +
  guides(colour=guide_legend(ncol=3, bycol =TRUE,title.position = 'left')) +
  xlab('Phycocyanin') + ylab(expression(paste(italic('p')*'CO'[2]~"(ppm)")))

## 3.  the variables which do not vary by day or night
## ===================================================================================================
targetlist <- c("airtemp", "stability", "dailyrain", "conv", "turb","cdom","windsp",
                "cond1")
realnames <- c("Chlorophyll", "Phycocyanin", "Air temperature", "Schmidt stability", 
               "Daily rain", "Convection", "Turbidity", "DOC",
               "Wind speed", "Conductivity (ca DIC)")
category <- c("Production", "Production", "Concentrations","Mixing","Mixing","Mixing","Decomposition",
              "Decomposition","Mixing","Concentrations")
varWantopt <- c("chlrfu", "bga1rfu", "airtemp", "stability", "dailyrain", "conv", "turb","cdom","windsp",
                "cond1") 

plotlist <- list()
for (i in 1: length(targetlist)) {
  target <- targetlist[[i]]
  N <- 200
  varWant <- varWantopt[- which(varWantopt == target)]
  lakeXbar <- as.data.frame(t(colMeans(bdat[, varWant], na.rm = TRUE)))
  
  take <- which(names(bdat) %in% target)
  co2.pdat <- data.frame(target = rep(seq(min(bdat[,take], na.rm = TRUE),
                                          max(bdat[,take], na.rm = TRUE),
                                          length = N), times=2),
                         TimeofDay = rep(c("Day", "Night"), each=N))
  co2.pdat$TimeofDay <- factor(co2.pdat$TimeofDay)
  names(co2.pdat)[1] <- target
  co2.pdat <- cbind(co2.pdat, lakeXbar)
  
  ## let's predict
  co2.pred <- predict(co2minus, newdata = co2.pdat, type = "terms", se.fit = TRUE)
  
  whichCols <- grep(target, colnames(co2.pred$fit))
  whichColsSE <- grep(target, colnames(co2.pred$se.fit))
  
  co2.resp <- cbind(co2.pdat, Fitted = co2.pred$fit[, whichCols], 
                    se.Fitted = co2.pred$se.fit[, whichColsSE])
  
  co2.resp <- with(co2.resp, transform(co2.resp, Fittedplus = Fitted + se.Fitted))
  co2.resp <- with(co2.resp, transform(co2.resp, Fittedminus = Fitted - se.Fitted))
  
  ## make into original limits
  shiftco2 <- attr(predict(co2minus, newdata = co2.pdat, type = "iterms"), "constant")
  co2.respnorm <- co2.resp
  co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fitted = Fitted + shiftco2))
  co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fittedplus = Fittedplus + shiftco2))
  co2.respnorm <- with(co2.respnorm, transform(co2.respnorm, Fittedminus = Fittedminus + shiftco2))
  whichexp <- grep("Fitted", names(co2.respnorm))
  co2.respnorm[,whichexp] <- exp(co2.respnorm[,whichexp])
  
  co2.respnorm <- co2.respnorm[-which(co2.respnorm$TimeofDay == "Night"),]
  co2.respnorm$target <- co2.respnorm[,which(names(co2.respnorm) == target)]
  
  bdat$target <- bdat[,which(names(bdat) == target)]
  realnameindex <- which(varWantopt %in% target)
  realname <- realnames[realnameindex]
  
  plotlist[[i]] <- ggplot(co2.respnorm, aes(x = target, y = Fitted)) +
    papertheme + 
    #annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
    geom_line() +
    geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                alpha = 0.25, fill='white', col="grey50") +  
    #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
    #         show.legend = FALSE, inherit.aes = FALSE) +
    #geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
    #geom_vline(xintercept = meanpH, linetype = 'dotted') +
    geom_rug3(aes(x=target, y=co2corr), data = bdat, stat = "identity", position = "identity", 
              sides = "bl", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
    theme(legend.position='top') +
    xlab(realname) + ylab(expression(paste(italic('p')*'CO'[2]~"(ppm)")))
}

#plotlist[c(1:length(plotlist))]

plot1 <- plotlist[[1]]
plot2 <- plotlist[[4]]
plot3 <- plotlist[[6]]
plot4 <- plotlist[[8]]
#plot_grid(cyanoplot, chlplot, plot1, plot2, plot3, plot4)


## by time plots
## Testing effect by time plots
## ==================================================================================
testing1 <- predict(co2minus, type = 'terms')
testing <- as.data.frame(testing1)
tosum <- grep("chl", colnames(testing))
chleffect <- rowSums(testing[,tosum], na.rm = TRUE)
testing <- testing[,-tosum]
testing$chlrfu <- chleffect

tosum <- grep("bga", colnames(testing))
cyanoeffect <- rowSums(testing[,tosum], na.rm = TRUE)
testing <- testing[,-tosum]
testing$bga1rfu <- cyanoeffect

names(testing)[1:8] <- varWantopt[-c(1:2)]

daynights <- testing[,c('chlrfu','bga1rfu')]
testing <- testing[,-grep('rfu', names(testing))]

## reorder bdat, since predictions all fucked up
testdf <- bdat[-which(bdat$bga1rfu >= 20),]
lasts <- testdf[which(!complete.cases(testdf[,c(varWantopt, 'co2corr')])),]
testdf <- testdf[which(complete.cases(testdf[,c(varWantopt, 'co2corr')])),]

testdfbyday <- testdf[order(testdf$TimeofDay, testdf$datetime),]
testdfbyday <- rbind(testdfbyday, lasts)

testdf <- rbind(testdf, lasts)

daynights$datetime <- testdfbyday$datetime
daynights$Year <- testdfbyday$Year
daynights$Hour <- testdfbyday$Hour
daynights$DOY <- testdfbyday$DOY
daynights$Month <- testdfbyday$Month
daynights$TimeofDay <- testdfbyday$TimeofDay

testing$datetime <- testdf$datetime
testing$Year <- testdf$Year
testing$Hour <- testdf$Hour
testing$DOY <- testdf$DOY
testing$Month <- testdf$Month
testing$TimeofDay <- testdf$TimeofDay


testsplit <- with(testing, split(testing, list(TimeofDay, DOY)))
testsplit <- testsplit[sapply(testsplit, function(x) dim(x)[1]) > 0] #remove empties for lakes R, E etc.

daynightsplit <- with(daynights, split(daynights, list(TimeofDay, DOY)))
daynightsplit <- daynightsplit[sapply(daynightsplit, function(x) dim(x)[1]) > 0] #remove empties for lakes R, E etc.

mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(TimeofDay = df['TimeofDay'][1,], DOY = df['DOY'][1,], t(means)))
}

testmeans <- do.call(rbind, lapply(testsplit, mycolMeans, cols=varWantopt[-c(1:2)]))
rownames(testmeans) <- NULL

daynightmeans <- do.call(rbind, lapply(daynightsplit, mycolMeans, cols=varWantopt[c(1:2)]))
rownames(daynightmeans) <- NULL

cyanos <- ggplot(daynightmeans, 
               aes(x=DOY, y=bga1rfu)) + #, fill=Lake, alpha=0.2
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.2, fill="grey60") +
  #geom_boxplot(show.legend = FALSE) +
  geom_line() +
  facet_wrap('TimeofDay', nrow=2) +
  ylab("Cyanobacteria effect") +
  theme(axis.title.x=element_blank())

chls <- ggplot(daynightmeans, 
               aes(x=DOY, y=chlrfu)) + #, fill=Lake, alpha=0.2
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.2, fill="grey60") +
  #geom_boxplot(show.legend = FALSE) +
  geom_line() +
  facet_wrap('TimeofDay', nrow=2) +
  ylab(expression(paste("Chl"~italic('a')~"effect"))) +
  theme(axis.title.x=element_blank())

conds <- ggplot(testmeans, 
               aes(x=DOY, y=cond1)) + #, fill=Lake, alpha=0.2
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.2, fill="grey60") +
  #geom_boxplot(show.legend = FALSE) +
  geom_line() +
  facet_wrap('TimeofDay', nrow=2) +
  ylab("Conductivity") +
  theme(axis.title.x=element_blank())

stabs <- ggplot(testmeans, 
               aes(x=DOY, y=stability)) + #, fill=Lake, alpha=0.2
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.2, fill="grey60") +
  #geom_boxplot(show.legend = FALSE) +
  geom_line() +
  facet_wrap('TimeofDay', nrow=2) +
  ylab("Stability effect") +
  theme(axis.title.x=element_blank())

convs <- ggplot(testmeans, 
               aes(x=DOY, y=conv)) + #, fill=Lake, alpha=0.2
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.2, fill="grey60") +
  #geom_boxplot(show.legend = FALSE) +
  geom_line() +
  facet_wrap('TimeofDay', nrow=2) +
  ylab("Potential convection effect") +
  theme(axis.title.x=element_blank())

cdoms <- ggplot(testing, 
               aes(x=DOY, y=cdom)) + #, fill=Lake, alpha=0.2
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.2, fill="grey60") +
  #geom_boxplot(show.legend = FALSE) +
  geom_line() +
  facet_wrap('TimeofDay', nrow=2) +
  ylab("Dissolved organic carbon effect") +
  theme(axis.title.x=element_blank())

airtemps <- ggplot(testmeans, 
               aes(x=DOY, y=airtemp)) + #, fill=Lake, alpha=0.2
  #annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.2, fill="grey60") +
  #geom_boxplot(show.legend = FALSE) +
  geom_line() +
  facet_wrap('TimeofDay', nrow=2) +
  ylab("Air temperature effect") +
  theme(axis.title.x=element_blank())


testmelt <- melt(testmeans, id.vars=c("TimeofDay", "DOY"))
testmelt <- merge(testmelt, data.frame(variable=varWantopt, realnames=realnames, category=category))

daynightmelt <- melt(daynightmeans, id.vars=c("TimeofDay", "DOY"))
daynightmelt <- merge(daynightmelt, data.frame(variable=varWantopt, realnames=realnames,
                                               category=category))

allmelt <- rbind(daynightmelt, testmelt)
allmelt$realnames <- factor(allmelt$realnames, levels = unique(allmelt$realnames)[c(3,5,4,9,6,7,8,10,
                                                                                    1,2)])

#saveRDS(allmelt, "../data/private/allmelt.rds")

ggplot(allmelt[-which(allmelt$variable == "chlrfu" | allmelt$variable == "windsp"),], 
       aes(y=value, x=DOY, col=category)) +
  papertheme +
  geom_line() +
  scale_color_viridis(name="Category", discrete=TRUE, end=0.9) +
  facet_grid(TimeofDay ~ realnames)

## predict for all vars in phmods
## FIXME: How can 2014 and 2015 be so different? e.g. conductivity non-overlap higher in 2015; many vars 
##    have completely opposite effect! cdom looks like in different units?? dailyrain crazy too; Conductivity
##    was higher in 2015 than any historically recorded conductivity by undergrads
## ============================================================================================================
## limit prediction data frame to observed intervals
regsplit <- with(bdat, split(bdat, list(TimeofDay)))
regsplit <- regsplit[sapply(regsplit, function(x) dim(x)[1]) > 0] #remove empties 
minmax <- function(df, colnames) {
  allmin <- as.data.frame(do.call(cbind, lapply(df[,colnames], min, na.rm = TRUE)))
  names(allmin) <- sapply(names(allmin), function(x) paste("min",x, sep = ""))
  allmax <- as.data.frame(do.call(cbind, lapply(df[,colnames], max, na.rm = TRUE)))
  names(allmax) <- sapply(names(allmax), function(x) paste("max",x, sep = ""))
  summ <- as.data.frame(cbind(allmin, allmax))
  summ <- cbind(summ, data.frame(TimeofDay = df$TimeofDay[1]))
  summ
}
minmaxes <- do.call(rbind, lapply(regsplit, minmax, colnames = c("chlrfu","bga1rfu", "airtemp", "stability", 
                                                                 "dailyrain", "conv", "turb","cdom",
                                                                 "windsp", "cond1")))
rownames(minmaxes) <- NULL

##  the variables which do not vary by day or night
## ===================================================================================================
targetlist <- c("airtemp", "stability", "dailyrain", "conv", "turb","cdom","windsp",
                "cond1")
realnames <- c("Chlorophyll", "Phycocyanin", "Air temperature", "Schmidt stability", 
               "Daily rain", "Convection", "Turbidity", "DOC",
               "Wind speed", "Conductivity (ca DIC)")
category <- c("Production", "Production", "Concentrations","Mixing","Mixing","Mixing","Decomposition",
              "Decomposition","Mixing","Concentrations")
varWantopt <- c("chlrfu", "bga1rfu", "airtemp", "stability", "dailyrain", "conv", "turb","cdom","windsp",
                "cond1") 

plotlist <- list()
for (i in 1: length(targetlist)) {
  target <- targetlist[[i]]
  N <- 200
  varWant <- varWantopt[- which(varWantopt == target)]
  lakeXbar <- as.data.frame(t(colMeans(bdat[, varWant], na.rm = TRUE)))
  
  take <- which(names(bdat) %in% target)
  ph.pdat <- data.frame(target = rep(seq(min(bdat[,take], na.rm = TRUE),
                                          max(bdat[,take], na.rm = TRUE),
                                          length = N), times=2),
                         TimeofDay = rep(c("Day", "Night"), each=N))
  ph.pdat$TimeofDay <- factor(ph.pdat$TimeofDay)
  names(ph.pdat)[1] <- target
  ph.pdat <- cbind(ph.pdat, lakeXbar)
  
  ## let's predict
  ph.pred <- predict(phmod, newdata = ph.pdat, type = "terms", se.fit = TRUE)
  
  whichCols <- grep(target, colnames(ph.pred$fit))
  whichColsSE <- grep(target, colnames(ph.pred$se.fit))
  
  ph.resp <- cbind(ph.pdat, Fitted = ph.pred$fit[, whichCols], 
                    se.Fitted = ph.pred$se.fit[, whichColsSE])
  
  ph.resp <- with(ph.resp, transform(ph.resp, Fittedplus = Fitted + se.Fitted))
  ph.resp <- with(ph.resp, transform(ph.resp, Fittedminus = Fitted - se.Fitted))
  
  ## make into original limits
  shiftco2 <- attr(predict(phmod, newdata = ph.pdat, type = "iterms"), "constant")
  ph.respnorm <- ph.resp
  ph.respnorm <- with(ph.respnorm, transform(ph.respnorm, Fitted = Fitted + shiftco2))
  ph.respnorm <- with(ph.respnorm, transform(ph.respnorm, Fittedplus = Fittedplus + shiftco2))
  ph.respnorm <- with(ph.respnorm, transform(ph.respnorm, Fittedminus = Fittedminus + shiftco2))
  whichexp <- grep("Fitted", names(ph.respnorm))
  
  ph.respnorm <- ph.respnorm[-which(ph.respnorm$TimeofDay == "Night"),]
  ph.respnorm$target <- ph.respnorm[,which(names(ph.respnorm) == target)]
  
  bdat$target <- bdat[,which(names(bdat) == target)]
  realnameindex <- which(varWantopt %in% target)
  realname <- realnames[realnameindex]
  
  plotlist[[i]] <- ggplot(ph.respnorm, aes(x = target, y = Fitted)) +
    papertheme + 
    #annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
    geom_line() +
    geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                alpha = 0.25, fill='white', col="grey50") +  
    #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
    #         show.legend = FALSE, inherit.aes = FALSE) +
    #geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
    #geom_vline(xintercept = meanpH, linetype = 'dotted') +
    geom_rug3(aes(x=target, y=ph1), data = bdat, stat = "identity", position = "identity", 
              sides = "bl", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
    theme(legend.position='top') +
    xlab(realname) + ylab("pH")
}


plotlist[c(1:length(plotlist))]

plotlist2015 <- list()
for (i in 1: length(targetlist)) {
  target <- targetlist[[i]]
  N <- 200
  varWant <- varWantopt[- which(varWantopt == target)]
  lakeXbar <- as.data.frame(t(colMeans(bdat5[, varWant], na.rm = TRUE)))
  
  take <- which(names(bdat5) %in% target)
  ph.pdat <- data.frame(target = rep(seq(min(bdat5[,take], na.rm = TRUE),
                                         max(bdat5[,take], na.rm = TRUE),
                                         length = N), times=2),
                        TimeofDay = rep(c("Day", "Night"), each=N))
  ph.pdat$TimeofDay <- factor(ph.pdat$TimeofDay)
  names(ph.pdat)[1] <- target
  ph.pdat <- cbind(ph.pdat, lakeXbar)
  
  ## let's predict
  ph.pred <- predict(phmod5, newdata = ph.pdat, type = "terms", se.fit = TRUE)
  
  whichCols <- grep(target, colnames(ph.pred$fit))
  whichColsSE <- grep(target, colnames(ph.pred$se.fit))
  
  ph.resp <- cbind(ph.pdat, Fitted = ph.pred$fit[, whichCols], 
                   se.Fitted = ph.pred$se.fit[, whichColsSE])
  
  ph.resp <- with(ph.resp, transform(ph.resp, Fittedplus = Fitted + se.Fitted))
  ph.resp <- with(ph.resp, transform(ph.resp, Fittedminus = Fitted - se.Fitted))
  
  ## make into original limits
  shiftco2 <- attr(predict(phmod, newdata = ph.pdat, type = "iterms"), "constant")
  ph.respnorm <- ph.resp
  ph.respnorm <- with(ph.respnorm, transform(ph.respnorm, Fitted = Fitted + shiftco2))
  ph.respnorm <- with(ph.respnorm, transform(ph.respnorm, Fittedplus = Fittedplus + shiftco2))
  ph.respnorm <- with(ph.respnorm, transform(ph.respnorm, Fittedminus = Fittedminus + shiftco2))
  whichexp <- grep("Fitted", names(ph.respnorm))
  
  ph.respnorm <- ph.respnorm[-which(ph.respnorm$TimeofDay == "Night"),]
  ph.respnorm$target <- ph.respnorm[,which(names(ph.respnorm) == target)]
  
  bdat5$target <- bdat5[,which(names(bdat5) == target)]
  realnameindex <- which(varWantopt %in% target)
  realname <- realnames[realnameindex]
  
  plotlist2015[[i]] <- ggplot(ph.respnorm, aes(x = target, y = Fitted)) +
    papertheme + 
    #annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
    geom_line() +
    geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                alpha = 0.25, fill='white', col="grey50") +  
    #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
    #         show.legend = FALSE, inherit.aes = FALSE) +
    #geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
    #geom_vline(xintercept = meanpH, linetype = 'dotted') +
    geom_rug3(aes(x=target, y=ph1), data = bdat5, stat = "identity", position = "identity", 
              sides = "bl", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
    theme(legend.position='top') +
    xlab(realname) + ylab("pH")
}


plotlist2015[c(1:length(plotlist2015))]

idx <- order(c(seq_along(plotlist), seq_along(plotlist2015)))
allplots <- c(plotlist,plotlist2015)[idx]

