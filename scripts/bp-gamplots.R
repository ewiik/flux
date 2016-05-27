## plots of the bp diel data

## load packages
library('ggplot2')
library('mgcv')

## read in models
co2mod <- readRDS('../data/private/bpco2phmod.rds')
phmod <- readRDS('../data/private/bpmod.rds')

## read in data
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- transform(bdat, daylength = sunset-sunrise)
bdat$daylength <- as.numeric(bdat$daylength)

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
  geom_vline(xintercept = meanpH, linetype="dotted") +
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
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanoxy, linetype="dotted") +
  xlab(expression(paste(O[2]~"%"))) + ylab('pH')

## save
allgam <- plot_grid(phplot, oxyplot, ncol = 2)

ggsave("../docs/private/bp-allgams.png", allgam, width=22, height=13, units = 'cm')
