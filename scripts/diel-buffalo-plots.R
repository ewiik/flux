## plots from modelling time on buffalo pound diel data

## load packages
library("mgcv")
library("ggplot2")
library("viridis")
library("extrafont")
library("gridExtra")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

## read in models
if(!file.exists("../data/private/BPphtimegammnull2015.rds")) {
  source("../scripts/diel-buffalo-models.R")}

gamm5null <- readRDS("../data/private/BPphtimegammnull2015.rds")
#gamm5AR1 <- readRDS("../data/private/BPphtimegammAR12015.rds")
## SFIXME: couldn't get the model to run in newest naming convention; too complex

#gamm4null <- readRDS("../data/private/BPphtimegammnull.rds") not created
gamm4AR1 <- readRDS("../data/private/BPphtimegammAR12014.rds")

bam5AR1 <- readRDS('../data/private/BPphtimebamAR12015.rds')


co2gammAR1 <- readRDS("../data/private/BPtimegammAR1.rds")
co2gammnull <- readRDS('../data/private/BPtimegammnull.rds')

## read in originals
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  source("diel-buffalo.R")}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- bdat[order(bdat$datetime),]

bdat5 <- readRDS('../data/private/bpbuoy2015-mod.rds')
bdat5 <- bdat5[order(bdat5$datetime),]

## 2014 & 2015 - let's check how different the pH is between the years!
##=======================================================================
master <- rbind(bdat[,c('ph1', 'ph2','DOY','datetime','Time', 'UpTime','DownTime')], 
                bdat5[,c('ph1','ph2','DOY','datetime', 'Time', 'UpTime','DownTime')])
master$Year <- format(master$datetime, format="%Y")

realplot <- ggplot(data=master, aes(y = Time, x = DOY, z=ph1, group=Year)) + 
  theme_bw() +
  facet_wrap('Year', nrow=2) +
  geom_raster(aes(fill=ph1)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent', name="pH") +
  geom_line(aes(y=UpTime)) +
  geom_line(aes(y=DownTime)) +
  ylab("Time of day; sunset & sunrise") +
  xlab("Day of year")

ggplot(data=master, aes(y = Time, x = DOY, z=ph1-ph2, group=Year)) + 
  theme_bw() +
  facet_wrap('Year', nrow=2) +
  geom_raster(aes(fill=ph1-ph2)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent') +
  geom_line(aes(y=UpTime)) +
  geom_line(aes(y=DownTime)) +
  ylab("Time of day; sunset & sunrise") +
  xlab("Day of year")
# biggest diffs when high ph....same controls?

pred2014 <- predict(gamm4AR1$gam, newdata=bdat, type = "response") # couldn't get NA rows to retain unless
#   specififed newdata...
pred2015 <- predict(gamm5AR1$gam, newdata=bdat5, type = "response")
pred2015null <- predict(gamm5null$gam, newdata=bdat5, type = "response")

allpred <- cbind(master, c(pred2014, pred2015))
names(allpred)[grep("pred",names(allpred))] <- "Fitted"
allpred$nullFitted <- rep(NA)
allpred$nullFitted[allpred$Year == 2015] <- pred2015null

Fittedplot <- ggplot(data=allpred, aes(y = Time, x = DOY, z=Fitted, group=Year)) + 
  theme_bw() +
  facet_wrap('Year', nrow=2) +
  geom_raster(aes(fill=Fitted)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent', name="Fitted pH") +
  geom_line(aes(y=UpTime)) +
  geom_line(aes(y=DownTime)) +
  ylab("Time of day; sunset & sunrise") +
  xlab("Day of year")

nullplot <- ggplot(data=allpred, aes(y = Time, x = DOY, z=nullFitted, group=Year)) + 
  theme_bw() +
  facet_wrap('Year', nrow=2) +
  geom_raster(aes(fill=nullFitted)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent', name="Fitted null pH") +
  geom_line(aes(y=UpTime)) +
  geom_line(aes(y=DownTime)) +
  ylab("Time of day; sunset & sunrise") +
  xlab("Day of year")

grid.arrange(realplot, Fittedplot, nullplot,ncol=3)
## FIXME: for 2015, the null model (no AR terms) is a lot better at summarising what happened...
##    The AR1 model just seems to miss the whole point..?!Why?

## 2014: 
## ============================================================================================
## quick basic models for lines for plots
lmnight <- with(bdat[bdat$isDay == FALSE,], lm(co2corr~pco2))
lmday <- with(bdat[bdat$isDay == TRUE,], lm(co2corr~pco2))
lmint <- with(bdat, lm(co2corr~pco2))
with(bdat, cor(x=co2corr, y=pco2, use='complete.obs', method='pearson'))

## the day time model plot
pglm <- ggplot(data=ten3, aes(y = co2corr, x = pco2)) + 
  theme_bw() +
  geom_point() +
  geom_abline(slope=1, intercept=0, color='red', lty=3) +
  #geom_ribbon() + ## FIXME: add the regression line here
  geom_abline(intercept = glmmod$coefficients[1], slope=glmmod$coefficients[2], col = 'black') +
  scale_color_viridis(discrete = TRUE, 
                      option = 'plasma', direction = -1) +
  ylab(expression(paste('Measured'~italic('p')*'CO'[2]~'('*mu*'atm)'))) +
  xlab(expression(paste('Calculated'~italic('p')*'CO'[2]~'('*mu*'atm)')))


## plot the difference between night and day in calculated vs measured.
with(bdat, plot(co2corr ~ pco2, col=ifelse(isDay, 'red', 'black'),
                xlab = 'Calculated pCO2', ylab = 'Measured pCO2'))
abline(0,1, lty=3)
abline(lmday$coefficients[1], lmday$coefficients[2], col = 'red')
abline(lmnight$coefficients[1], lmnight$coefficients[2], col = 'black')
#abline(-71.6, 1.143408, lty = 3, col = 'green')
legend('topleft', col = c('red', 'black', 'black'), pch=c(1,1, NA), lty = c(1,1,3), 
       legend = c('daytime values and regression line', 
                  'night-time values and regression line',
                  '1:1 line'))
legend('bottomright', col = c('red', 'black'), title = 'R2', 
       legend = c('day = 0.72', 'night = 0.92'))

with(bdat, plot(pco2 - co2corr ~ co2corr, xlab = 'Measured pCO2', 
                ylab = 'Calculated - measured pCO2'))
abline(0,0)


##predict for periods of varying interactions.
N <- 200
simDOY <- c(171:173, 198:200, 237:239)
DOYgroup <- factor(rep(c('Jun 20-22', 'Jul 17-19', 'Aug 25-27'), times=c(3,3,3)))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(bdat$Time, na.rm=TRUE),
                                      max(bdat$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
mpred <- predict(gamm4AR1$gam, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, mpred, DOYgroups)
predicted$DOYgroups <- factor(predicted$DOYgroups, levels = unique(predicted$DOYgroups)[c(1,2,3)])
names(predicted)[which(names(predicted)=='mpred')] <- 'pH'

predph2014 <- predicted
ph2014 <- ggplot(predicted, aes(x = Time, y = pH, group= DOY, col=DOYgroups)) +
  papertheme +
  geom_line() +
  scale_colour_brewer(name = "Date", type = 'seq', palette = 'PuBuGn', direction=1) +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab('pH')

## for co2
N <- 200
simDOY <- c(168:170, 190:192, 228:230)
DOYgroup <- factor(rep(1:3, times=c(3,3,3)))
DOYgroup <- factor(rep(c('Jun 17-19', 'Jul 9-11', 'Aug 16-18'), times=c(3,3,3)))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(bdat2014full$Time, na.rm=TRUE),
                                      max(bdat2014full$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
mco2pred <- predict(co2gammAR1$gam, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, mco2pred, DOYgroups)
names(predicted)[which(names(predicted)=='mco2pred')] <- 'pCO2'

labdatco2 <- data.frame(x = 3, y = 380, label = "Mean (atm)")

co2modplot <- ggplot(predicted, aes(x = Time, y = pCO2, group= DOY, col=DOYgroups)) +
  theme_bw(base_size = 10) +
  scale_colour_brewer(name = "Date", type = 'qual', palette = 'Dark2', direction=1) +
  geom_line() +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab(expression(paste(italic(p)*"CO"[2]~"(ppm)"))) +
  geom_abline(intercept = 398, slope = 0, linetype='dotted') +
  geom_text(data = labdatco2, aes(label = label, x = x, y = y), 
            show.legend = FALSE, inherit.aes = FALSE,  size = 3)
ggsave('../docs/private/bp-diel-co2gam.png', width=20, height=15, units='cm')

## arrange a few plots
asloplots <- plot_grid(co2modplot, dielplot, ncol=2,label_size = 10)
ggsave('../docs/private/dielplots.png', asloplots, width=20, height=10, units='cm')

## predict for the gamm models 
## (as per http://www.fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/)
want <- seq(1, nrow(bdat), length.out = 200)
pdat <- with(bdat, data.frame(Time = Time[want], DOY = DOY[want]))
 
## predict trend contributions
p  <- predict(gamm4null$gam,  newdata = pdat, type = "terms", se.fit = TRUE)
p1 <- predict(gamm4AR1$gam, newdata = pdat, type = "terms", se.fit = TRUE)

## combine with the predictions data, including fitted and SEs
pdat <- transform(pdat, p  = p$fit[,1],  se  = p$se.fit[,1],
                  p1 = p1$fit[,1], se1 = p1$se.fit[,1])

op <- par(mar = c(5,4,2,2) + 0.1)
ylim <- with(pdat, range(p, p1))
ylim[1] <- floor(ylim[1])
ylim[2] <- ceiling(ylim[2])
ylab <- "pCO2"
plot(co2corr - mean(co2corr, na.rm = TRUE) ~ DOY, data = bdat, type="l", lty=2,col="green",
      ylab = ylab, ylim = ylim)
lines(p  ~ DOY, data = pdat, col = "black")
points(p1 ~ DOY, data = pdat, col = "red")
lines(co2corr - mean(co2corr) ~ DOY, data = bdat, lty=2,col="green")
legend("topleft",
        legend = c("Uncorrelated Errors", paste0("AR(", 1, ") Errors")),
        bty = "n", col = c("black","red"),
        lty = 1, lwd = c(1,1,1))
par(op)

## ============================================================================================
## for 2015
## ============================================================================================
N <- 200
simDOY <- c(159:161, 190:192, 237:239)
DOYgroup <- factor(rep(c('Jun 8-10', 'Jul 9-11', 'Aug 25-27'), times=c(3,3,3)))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(bdat5$Time, na.rm=TRUE),
                                      max(bdat5$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
mpred <- predict(gamm5AR1$gam, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, mpred, DOYgroups)
names(predicted)[which(names(predicted)=='mpred')] <- 'pH'
predicted$DOYgroups <- factor(predicted$DOYgroups, levels = unique(predicted$DOYgroups)[c(1,2,3)])
predph2015 <- predicted

ph2015 <- ggplot(predicted, aes(x = Time, y = pH, group= DOY, col=DOYgroups)) +
  papertheme +
  geom_line() +
  #scale_colour_brewer(name = "Date", type = 'seq', palette = 'PuBuGn', direction=1) +
  scale_color_viridis(name = "Date", discrete = TRUE, begin=0.2, end=0.8) +
  geom_line() +
  xlab('Time of Day') + ylab('pH')

## for 2014 and 2015 together
predph2015$Year <- rep(2015)
predph2014$Year <- rep(2014)
predph <- rbind(predph2014, predph2015)
predph$DOYgroups <- factor(predph$DOYgroups, levels = unique(predph$DOYgroups)[c(4,1,5,2,3)])

ggplot(predph, aes(x = Time, y = pH, group= DOY, col=DOYgroups)) +
  papertheme +
  geom_line() +
  #scale_colour_brewer(name = "Date", type = 'seq', palette = 'PuBuGn', direction=1) +
  scale_color_viridis(name = "Date", option="inferno", discrete = TRUE, begin=0.2, end=0.8) +
  facet_wrap('Year') +
  geom_line() +
  xlab('Time of Day') + ylab('pH')