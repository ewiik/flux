## plots of the bp diel data

## load packages
library('ggplot2')
library('mgcv')
library('cowplot')

## read in models
co2mod <- readRDS('../data/private/bpco2phmod.rds')
phmod <- readRDS('../data/private/bpmod.rds')

bptimemod <- readRDS("../data/private/bp-diel-timemod2014.rds")
bpco2mod <- readRDS("../data/private/bpmod.rds")

## read in data
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- transform(bdat, daylength = sunset-sunrise)
bdat$daylength <- as.numeric(bdat$daylength)

## predict for all vars in bpco2mod
## 1. chl 
N <- 200
varWant <- c("bga1rfu", "airtemp", "stability", "ODOrel1", "dailyrain", "conv", "turb","cdom","windsp")
bdat$TimeofDay <- as.factor(bdat$TimeofDay)
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
co2.pred <- predict(bpmod, newdata = co2.pdat, type = "terms", se.fit = TRUE)
co2.test <- predict(bpmod, newdata = co2.pdat, type = "response", se.fit = TRUE)
whichCols <- grep("chl", colnames(co2.pred$fit))
whichColsSE <- grep("chl", colnames(co2.pred$se.fit))
co2.pdat <- cbind(co2.pdat, Fitted = rowSums(co2.pred$fit[, whichCols]), 
                  se.Fitted = rowSums(co2.pred$se.fit[, whichColsSE]))

co2.pdat <- with(co2.pdat, transform(co2.pdat, Fittedplus = Fitted + se.Fitted))
co2.pdat <- with(co2.pdat, transform(co2.pdat, Fittedminus = Fitted - se.Fitted))

## make into original limits
shiftco2 <- attr(predict(bpmod, newdata = co2.pdat, type = "iterms"), "constant")
co2.pdatnorm <- co2.pdat
co2.pdatnorm <- with(co2.pdatnorm, transform(co2.pdatnorm, Fitted = Fitted + shiftco2))
co2.pdatnorm <- merge(co2.pdatnorm, minmaxes)
overs <- with(co2.pdatnorm, which(pH_surface < minpH_surface | pH_surface > maxpH_surface))
co2.pdatnorm <- co2.pdatnorm[-overs,]

labdatco2 <- data.frame(x = 7.5, y = meanco2 - 18, label = "mean flux")

## add quantiles
phquants <- quantile(regvarf2$pH_surface, c(.05,.95), na.rm = TRUE)

co2plot <- ggplot(co2.pdatnorm, aes(x = pH_surface, y = Fitted, 
                                    colour = ifelse(Lake == "WW", "Wascana", 
                                                    ifelse(Lake == "D", "Diefenbaker",
                                                           ifelse(Lake == 'K', "Katepwa", 
                                                                  ifelse(Lake == 'P', "Pasqua", 
                                                                         ifelse(Lake == 'B', 'Buffalo Pound', 
                                                                                ifelse(Lake=='L','Last Mountain',
                                                                                       'Crooked')))))),
                                    lty=ifelse(Lake == "WW", "Wascana", 
                                               ifelse(Lake == "D", "Diefenbaker",
                                                      ifelse(Lake == 'K', "Katepwa", 
                                                             ifelse(Lake == 'P', "Pasqua", 
                                                                    ifelse(Lake == 'B', 'Buffalo Pound', 
                                                                           ifelse(Lake=='L','Last Mountain',
                                                                                  'Crooked')))))))) +
  papertheme + 
  annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  #geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
  #           alpha = 0.25, fill='white') +  
  scale_linetype_manual(name='Lake', values = c("solid", "dotdash","longdash", "solid", "longdash", 
                                                "solid", "solid")) +
  scale_colour_manual(name="Lake", values = c("#b2abd2", "#5e3c99","#b2abd2", "#5e3c99", "#e66101",
                                              "#5e3c99", "#5e3c99"))+
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE, inherit.aes = FALSE) +
  geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  geom_vline(xintercept = meanpH, linetype = 'dotted') +
  geom_rug3(aes(x=pH_surface, y=co2Flux), data = regvarf2, stat = "identity", position = "identity", 
            sides = "bl", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  theme(legend.position='top') +
  guides(colour=guide_legend(ncol=3, bycol =TRUE,title.position = 'left')) +
  xlab('pH') + ylab(expression(paste('CO'[2]*' (mmol m'^{-2}*d^{-1}*')')))



## ==============================================================================================
## everything below this is unchecked
## ==============================================================================================
## means
meanco2 <- mean(bdat$co2corr, na.rm=TRUE)
meanph <- mean(bdat$ph1, na.rm = TRUE)
meanoxy <- mean(bdat$ODOrel1, na.rm = TRUE)

## predict for co2mod
N <- 200
null.pdat <- with(bdat,
                  data.frame(`ph1` = seq(min(ph1, na.rm=TRUE),
                                             max(ph1, na.rm=TRUE),length = N)))
null.pred <- predict(co2mod, newdata = null.pdat, type = "iterms") 
whichCols <- grep("ph", colnames(null.pred))

null.predse <- predict(co2mod, newdata = null.pdat, type = "iterms", se.fit = TRUE)
null.predse <- as.data.frame(null.predse$se.fit)
whichColsse <- grep("ph", colnames(null.predse))

null.pdat <- cbind(null.pdat, Fitted = null.pred[, whichCols], Fittedse = null.predse[,whichColsse])
null.pdat <- with(null.pdat, transform(null.pdat, Fittedplus = Fitted + Fittedse))
null.pdat <- with(null.pdat, transform(null.pdat, Fittedminus = Fitted - Fittedse))

shiftnull <- attr(predict(co2mod, newdata = null.pdat, type = "iterms"), "constant")
null.pdatnorm <- null.pdat
null.pdatnorm <- with(null.pdatnorm, transform(null.pdatnorm, Fitted = Fitted + shiftnull, 
                                               Fittedplus = Fittedplus + shiftnull, 
                                               Fittedminus = Fittedminus + shiftnull))

labdatco2 <- data.frame(x = 7.5, y = meanco2 - 18, label = "mean flux")

phplot <- ggplot(null.pdatnorm, aes(x = ph1, y = Fitted)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
  #  show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  geom_vline(xintercept = meanph, linetype="dotted") +
  ylab(expression(paste(italic(p)*"CO"[2]~"(ppm)"))) + xlab('pH')

## predict for phmod
N <- 200
varWant <- c("chl", "airtemp", "daylength")
bdatsub <- bdat[,varWant]
lakeXbar <- colMeans(bdatsub, na.rm = TRUE)
lakeXbar <- t(as.data.frame(lakeXbar))
rownames(lakeXbar) <- NULL

oxy.pdat <- with(bdat,
                 data.frame(`ODOrel1` = seq(min(`ODOrel1`, na.rm = TRUE),
                                                   max(`ODOrel1`, na.rm = TRUE),
                                                   length = N)))
oxy.pdat <- cbind(oxy.pdat, lakeXbar)
oxy.pred <- predict(phmod, newdata = oxy.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("ODO", colnames(oxy.pred$fit))
whichColsSE <- grep("ODO", colnames(oxy.pred$se.fit))
oxy.pdat <- cbind(oxy.pdat, Fitted = oxy.pred$fit[, whichCols], 
                  se.Fitted = oxy.pred$se.fit[, whichColsSE])
## make into original limits
oxy.pdat <- with(oxy.pdat, transform(oxy.pdat, Fittedplus = Fitted + se.Fitted))
oxy.pdat <- with(oxy.pdat, transform(oxy.pdat, Fittedminus = Fitted - se.Fitted))

shiftoxy <- attr(predict(phmod, newdata = oxy.pdat, type = "iterms"), "constant")
oxy.pdatnorm <- oxy.pdat
oxy.pdatnorm <- with(oxy.pdatnorm, transform(oxy.pdatnorm, Fitted = Fitted + shiftoxy, 
                                             Fittedplus = Fittedplus + shiftoxy, 
                                             Fittedminus = Fittedminus + shiftoxy))
labdatoxy <- data.frame(x = c(2.5, 11), y = c(meanpH + 0.03, 8.4), label = c("mean pH", 'mean oxygen'))

oxyplot <- ggplot(oxy.pdatnorm, aes(x = ODOrel1, y = Fitted)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanph, linetype="dotted") +
  geom_vline(xintercept = meanoxy, linetype="dotted") +
  xlab(expression(paste(O[2]~"%"))) + ylab('pH')

## save
allgam <- plot_grid(phplot, oxyplot, ncol = 2)

ggsave("../docs/private/bp-allgams.png", allgam, width=22, height=13, units = 'cm')
