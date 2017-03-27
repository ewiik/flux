## statistics on bp data
stop("This script not ready to be blindly sourced yet! Execute chunkwise!")

## Notes on measures and instrumentation: CDOM is based on fluorometry:
##    http://www.turnerdesigns.com/t2/doc/appnotes/S-0022.pdf
##    Turbidity is this: YSI 6136 Turbidity Sensor; see this http://or.water.usgs.gov/grapher/fnu.html
## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
    source("diel-buffalo.R")         # creates both 2014 and 2015 data
}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- bdat[order(bdat$datetime),]

bdat5 <- readRDS('../data/private/bpbuoy2015-mod.rds')
bdat5 <- bdat5[order(bdat5$datetime),]

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

## load packages
library("mgcv")
library("ggplot2")
library("viridis")

## create ARn term function
myAR <- function(term, lag) {
  lagterm <- vector(length=length(term))
  for(i in 1:length(term)) {
    if(i <= lag) { lagterm[i] <- NA} else {
      lagterm[i] <- term[i-lag]
    }}
  lagterm
}

## ==================================================================================
## 2014: all data available
## ==================================================================================

## are there differences in relationship between night and day?
## ===================================================================================
## simple correlations -- ok for correlation coefficients but not for predictions without
##    taking care of heteroscedasticity
lmnight <- with(bdat[bdat$isDay == FALSE,], lm(co2corr~pco2))
lmday <- with(bdat[bdat$isDay == TRUE,], lm(co2corr~pco2))
lmint <- with(bdat, lm(co2corr~pco2))
with(bdat, cor(x=co2corr, y=pco2, use='complete.obs', method='pearson'))

## just with daytime data
days <- subset(bdat, isDay == TRUE)
keep <- which(complete.cases(days[,c('pco2','co2corr')]))

## try a glm and look at residuals etc.
daynull <- glm(co2corr ~ pco2, data = days[keep,], family = gaussian, x=TRUE)
plot(resid(daynull) ~ daynull$fitted.values, ylab='Residuals (pCO2 uatm)',
     xlab='Fitted values (pCO2 uatm)')
plot(resid(daynull) ~ daynull$y, ylab='Residuals (pCO2 uatm)',
     xlab='Response values (pCO2 uatm)')
plot(resid(daynull) ~ daynull$x[,'pco2'], ylab='Residuals (pCO2 uatm)',
     xlab='Predictor values (pCO2 uatm)')
## will boxcoxing help?
boxed <- boxcox(daynull)
daybox <- glm(co2corr^0.47 ~ pco2, data = days[keep,], family = gaussian, x=TRUE) #no

## will gls help?
## FIXME: Gavin, all the lit on gls and variance functions were really densely and unhelpfully
##    written --- haven't gotten further than this - no idea really how to test the best varFunc
vf1 <- varPower(form=~pco2)
vf1 <- Initialize(vf1, bdat)
coef(vf1) <- 0.3
varWeights(vf1)[1:10]

dayvar <- gls(co2corr ~ pco2, data = days[keep,], weights = vf1)
plot(dayvar)

## gamming the time component of the data for CO2
## FIXME: pH models still undeveloped
## =================================================================================================
## create lag terms
bdat <- transform(bdat, timelag1 = myAR(bdat$Time, 1))
bdat <- transform(bdat, co2lag1 = myAR(bdat$co2corr, 1))

## create models; testing AR inclusion and comparing k with edf to set k (gam.check)
## note that "in recent mgcv versions the output in summary() is more reliable than
##    the generalised likelihood ratio test we might do to compare a selected model
##    with a null one" (GS Oct 16)
## FIXME: DOY is eating up as many k's as possible; 80 too many, using most of 70
## FIXME: update model by Gav's ryver comments!
mco2ti <- gam(co2corr ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=70) + s(co2lag1) +
              ti(Time, DOY, bs = c("cc","tp")), # worried of setting ti k higher, kept crashing
            data = bdat, select = TRUE, method = "REML", family = scat,
            na.action = na.exclude)
# interaction term significant -->
mco2te <- gam(co2corr ~ te(Time, DOY, bs = c("cc","tp")) + s(co2lag1),
              data = bdat, select = TRUE, method = "REML", family = scat,
            na.action = na.exclude)
## FIXME: changing k makes it crash??? Does te() handle k differently?
## co2 dependence on previous value dwarfs all other terms... not sure I'm including everything
##    correctly? was i meant to do the autoregression separately?

plot(acf(resid(mco2te), na.action = na.exclude)) # still autocorrelated.. could add AR2?

## what if we add Week as a variable too? How to resolve interactions in this case?
mco2full <- gam(co2corr ~ te(Time, DOY, bs = c("cc", "tp")) + s(Week), data=bdat,
                select = TRUE,
                method = "REML", family = scat, na.action = na.exclude)

## try gamm approach following fromthebottomoftheheap?
## FIXME: should we pursue this?
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
testingnull <- gamm(co2corr ~te(Time, DOY, bs = c("cc","tp")), data = bdat,
                    control = ctrl, verbosePQL = TRUE)
testing1 <- gamm(co2corr ~te(Time, DOY, bs = c("cc","tp")), data = bdat, 
                 correlation = corARMA(form = ~ 1, p=1),
                 control = ctrl, verbosePQL = TRUE)
##testing2 <- gamm(co2corr ~te(Time, DOY, bs = c("cc","tp")), data = bdat, 
       ##          correlation = corARMA(form = ~ 1, p=2),
      ##           control = ctrl, verbosePQL = TRUE)
## FIXME: "Error in `coef<-.corARMA`(`*tmp*`, value = value[parMap[, i]]) : 
##    Coefficient matrix not invertible""
res <- resid(testing1$lme, type = "normalized")
acf(res, lag.max = 200, main = "ACF - AR(2) errors")

phgamm <- gamm(ph1 ~te(Time, DOY, bs = c("cc","tp")), data = bdat, 
                 correlation = corARMA(form = ~ 1, p=1),
                 control = ctrl, verbosePQL = TRUE)


## gam on other params: based on discussions with Helen, we should use rfu not mg/L
## ===========================================================================================
mbio <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay)) + s(bga1rfu, by = factor(TimeofDay)) + 
              s(airtemp) + s(stability) + s(ODOrel1) + s(dailyrain) + s(conv) + 
              s(turb) + s(cdom) + s(windsp) + s(cond1),
            data = bdat[bdat$bga1rfu < 20,], select = TRUE, method = "REML", family = tw,
            na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
## FIXME - see below , should oxygen also be removed?
## FIXME - what is up with bga1rfu? skews predictions massivley

mbiominus <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay), k=4) + 
                   s(bga1rfu, by = factor(TimeofDay), k=4) + 
              s(airtemp) + s(log(stability+1)) + s(log(dailyrain+1), k=4) + s(conv) + 
              s(log(turb), k=4) + s(cdom) + s(windsp) + s(cond1, k=4),
            data = bdat[bdat$bga1rfu < 20,], select = TRUE, method = "REML", family = tw,
            na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))

mbiominustest <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay), k=4) + 
                   s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
                   s(airtemp) + s(log(stability+1), k=3) + s(log(dailyrain+1), k=4) + s(conv) + 
                   s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
                 data = bdat, select = TRUE, method = "REML", family = tw,
                 na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
# had to put k for some vars to reduce overfitting
## no oxygen/ph here since too much sameness going on between co2, o2 and ph

mbiominusti <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay), k=4) + 
                     s(bga1rfu, by = factor(TimeofDay), k=4) + 
                     s(airtemp) + ti(log(stability+1)) + s(log(dailyrain+1), k=4) + ti(conv) +
                     s(log(turb), k=4) + s(cdom) + s(windsp) + s(cond1, k=4),
                   data = bdat[bdat$bga1rfu < 20,], select = TRUE, method = "REML", family = tw,
                   na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))

mbiominuste <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay), k=4) + 
                   s(bga1rfu, by = factor(TimeofDay), k=4) + 
                   s(airtemp) + ti(log(stability+1)) + s(log(dailyrain+1), k=4) + ti(conv) +
                   ti(log(stability +1), conv) +
                   s(log(turb), k=4) + s(cdom) + s(windsp) + s(cond1, k=4),
                 data = bdat[bdat$bga1rfu < 20,], select = TRUE, method = "REML", family = tw,
                 na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
anova(mbiominusti, mbiominuste, test="LRT")
## according to anova the interaction isn't important..
## FIXME -- gam's p-values have seemed a bit suspicious for a while... higher signif than expected
##    for many things... has something changed in computation? (e.g. windsp, with v low F)

## what about modeling pH since in fact these other vars might actually explain both?
phbio <- gam(ph1 ~ s(log(chlrfu), by = factor(TimeofDay), k=4) + s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
              s(airtemp) + s(log(stability +1), k=3) + s(log(dailyrain+1), k=4) + s(conv) + 
              s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
            data = bdat, select = TRUE, method = "REML", family = gaussian(),
            na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
# had to put k for some vars to reduce overfitting

## modeling co2 only based on oxygen and pH separated again since explains most everything
co2phmod <- gam(co2corr ~ s(ph1) + s(ODOrel1), method = "REML", family = tw(),
                na.action = na.exclude, data=bdat, select=TRUE,
                control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
co2gamma <- gam(co2corr ~ ti(ph1) + ti(ODOrel1) + ti(ph1, ODOrel1), method = "REML", family = Gamma(),
                na.action = na.exclude, data=bdat,
                control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
co2twee<- gam(co2corr ~ ti(ph1) + ti(ODOrel1) + ti(ph1, ODOrel1), method = "REML", family = tw(),
              na.action = na.exclude, data=bdat,
              control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
## FIXME -- why is the interaction significant when the relationship really just looks simple?
## FIXME -- is tweedie the best one here?

## save models/output
## =========================================================================================
pdf("../docs/private/diel-bpgam.pdf")
opar <- par()
par(mar = c(4,4,1,1) + 0.1)
plot(mco2te, pages = 1, scheme = 2)
par(opar)
dev.off()

saveRDS(bpmod, '../data/private/bpmod.rds')
saveRDS(mbiominustest, '../data/private/bpmodsimp.rds')
saveRDS(bpco2phmod, '../data/private/bpco2phmod.rds')
saveRDS(testingnull, '../data/private/BPgammnull.rds')
saveRDS(testing1, '../data/private/BPgammAR1.rds')
saveRDS(phbio, '../data/private/BPphmod.rds')
saveRDS(phgamm, '../data/private/BPphgamm.rds')

## ============================================================================================
## GAVIN MAY STOP READING HERE!
## ============================================================================================

## ==========================================================================================
## 2015: no CO2 or O2 shallow data but other params ok
## ==========================================================================================
## are the controls of pH the same in 2015?
phbio5 <- gam(ph1 ~ s(log(chlrfu), by = factor(TimeofDay), k=4) + s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
               s(airtemp) + s(log(stability +1), k=3) + s(log(dailyrain+1), k=4) + s(conv) + 
               s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
             data = bdat5, select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))

## gamming the time component of the data
m1 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(DOY),
          data = bdat5)
m2 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")),
          data = bdat5)
anova(m1, m2, test = "LRT") # m2 better but not convinced it's significant?
m3 <- gam(ph1 ~ te(Time, DOY, bs=c("cc", "tp")), data=bdat5)

mco2 <- gam(co2corr ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")),
            data = bdat5)

plot(m2, scheme=2)
plot(acf(resid(m3)))

## try gamm approach following fromthebottomoftheheap?
## FIXME: should we pursue this?
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
testingnull5 <- gamm(ph1 ~te(Time, DOY, bs = c("cc","tp")), data = bdat5,
                    control = ctrl, verbosePQL = TRUE)
testing15 <- gamm(ph1 ~te(Time, DOY, bs = c("cc","tp")), data = bdat5, 
                 correlation = corARMA(form = ~ 1, p=1),
                 control = ctrl, verbosePQL = TRUE) #for this model to run, I had to kill firefox
#   "Error in recalc.corAR1(cSt, list(Xy = as.matrix(val))) : 
#   'Calloc' could not allocate memory (305795169 of 8 bytes)""

##testing2 <- gamm(co2corr ~te(Time, DOY, bs = c("cc","tp")), data = bdat, 
##          correlation = corARMA(form = ~ 1, p=2),
##           control = ctrl, verbosePQL = TRUE)
## FIXME: Haven't tried this yet since for 2014 data it doesn't fit with matrix error

res <- resid(testing15$lme, type = "normalized")
resnull <- resid(testingnull5$lme, type = "normalized")

op <-par(mfrow=c(2,1))
acf(resnull, lag.max = 200, main = "ACF - default errors")
acf(res, lag.max = 200, main = "ACF - AR(1) errors")
par(op)


## save best models
saveRDS(phbio5, "../data/private/BPphmod2015.rds")
saveRDS(testing15, "../data/private/bp-diel-timemod2015.rds")
saveRDS(testingnull5, "../data/private/bp-diel-nullmod2015.rds")
