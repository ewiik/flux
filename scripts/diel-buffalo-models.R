## statistics on bp data
stop("This script not ready to be blindly sourced yet! Execute chunkwise!")

## Notes on measures and instrumentation: CDOM is based on fluorometry: 
##    http://www.turnerdesigns.com/t2/doc/appnotes/S-0022.pdf
##    Turbidity is this: YSI 6136 Turbidity Sensor; see this http://or.water.usgs.gov/grapher/fnu.html
## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  source("diel-buffalo.R")} #creates both 2014 and 2015 data
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- bdat[order(bdat$datetime),]

bdat5 <- readRDS('../data/private/bpbuoy2015-mod.rds')
bdat5 <- bdat5[order(bdat5$datetime),]

if (!file.exists("../data/private/bp-stability.rds")) {
  source("../scripts/bp-stratification.R")
}
schmidt <- readRDS("../data/private/bp-stability.rds") # note this for 2014
bdat <- merge(bdat, schmidt)

## create index of potential convective mixing, based on diff in T btw air and water
bdat$conv <- bdat$airtemp - bdat$temp4 # this is the topmost, 77cm one

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
testing2 <- gamm(co2corr ~ s(Time, bs = "cc", k = 20) + s(DOY, k = 30),
                 data = bdat, correlation = corCAR1(form = ~ datetime), 
                 control = ctrl)
## gamm crashes if I specify 0.8 as the correlation! in fact gamm crashes almost with anything
##    I try!
res <- resid(testing2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")

## for pH... To be developed! (ignore section)
##==========================================================================================
mnull <- gam(ph1 ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=30), data=bdat) 
m1 <- gam(ph1 ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=30) + ti(Time, DOY, bs = c("cc","tp")), 
          data = bdat)
m2 <- gam(ph1 ~ te(Time, DOY, bs = c("cc","tp")), data=bdat)

plot(acf(resid(m1)))

anova(mnull, m1, test = "LRT") # m1 better

## gam on other params: based on discussions with Helen, we should use rfu not mg/L
## ===========================================================================================
bpmod <- gam(co2corr ~
               s(chlrfu) + 
               s(bga1rfu) +
               s(airtemp) +
               s(stability) +  
               s(ODOrel1) +
               s(ph1)  +
               s(dailyrain) +
               s(conv) +
               s(turb) +
               s(windsp), 
             data = bdat,
             select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE,
                                   newton = list(maxHalf = 60)))
## FIXME: tested scat but though histogram looked better nothing else did...??
## FIXME: O2 nuked by pH. What to do about this? (tested interaction visually as you suggested
##    and no O2-pH-CO2 pattern seemed apparent)

## pH again can explain most everything
bpco2phmod <- gam(co2corr ~ s(ph1), method = "REML", family = gaussian,
                  na.action = na.exclude, data=bdat,
                  control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))

## what about modeling pH since in fact these other vars might actually explain both?
phmod <- gam(ph1 ~
               s(chlrfu) + 
               s(bga1rfu) +
               s(airtemp) +
               s(stability) +  
               s(ODOrel1) +
               s(dailyrain) +
               s(conv) +
               s(turb) +
               s(windsp), 
             data = bdat,
             select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE,
                                   newton = list(maxHalf = 60)))

## save models/output
## =========================================================================================
pdf("../docs/private/diel-bpgam.pdf")
opar <- par()
par(mar = c(4,4,1,1) + 0.1)
plot(mco2te, pages = 1, scheme = 2)
par(opar)
dev.off()

saveRDS(bpmod, '../data/private/bpmod.rds')

saveRDS(bpco2phmod, '../data/private/bpco2phmod.rds')

## ============================================================================================
## GAVIN MAY STOP READING HERE!
## ============================================================================================

## ==========================================================================================
## 2015: no CO2 or O2 shallow data but other params ok
## ==========================================================================================

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

saveRDS(m3, "../data/private/bp-diel-timemod2015.rds")
