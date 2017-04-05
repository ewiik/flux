## Time and mechanistic models on Buffalo Pound diel data (2014 and 2015)
## First, time and mechanistic models of CO2 in 2014
## Second, same models but for pH, for 2014 and 2015

## note on gamma models: http://stats.stackexchange.com/questions/67547/when-to-use-gamma-glms

## Notes on measures and instrumentation: CDOM is based on fluorometry:
##    http://www.turnerdesigns.com/t2/doc/appnotes/S-0022.pdf
##    Turbidity is this: YSI 6136 Turbidity Sensor; see this http://or.water.usgs.gov/grapher/fnu.html

## load packages
library("mgcv")
library("ggplot2")
library("viridis")
library("itsadug") # residual visualisation for bams with AR 
library("MASS") # for boxcox function

## are we running absolutely all models or just final models?
runextras <- FALSE

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
bdat5$conv <- bdat5$airtemp - bdat5$temp4 # this is the 77cm one (temp4 at 45cm)


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

## Measured vs calculated pCO2: are there differences in relationship between night and day?
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

if(runextras) {
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
}

## pCO2 vs time (2014)
## =================================================================================================
if(runextras){
  ## create lag terms
  bdat <- transform(bdat, timelag1 = myAR(bdat$Time, 1))
  bdat <- transform(bdat, co2lag1 = myAR(bdat$co2corr, 1))
  
  ## create models; testing AR inclusion and comparing k with edf to set k (gam.check)
  ## note that "in recent mgcv versions the output in summary() is more reliable than
  ##    the generalised likelihood ratio test we might do to compare a selected model
  ##    with a null one" (GS Oct 16)
  ## FIXME: DOY is eating up as many k's as possible; 80 too many, using most of 70
  ## FIXME: update model by Gav's ryver comments! "Made some progress with time model. 
  ## Change "tp" to "gp" for both DoY terms and remove any k values. Use Gamma and log link. 
  ##    Model will be ok-ish for now. Will look at model with predictors in AM"
  mco2ti <- gam(co2corr ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=70) + s(co2lag1) +
                  ti(Time, DOY, bs = c("cc","tp")), # worried of setting ti k higher, kept crashing
                data = bdat, select = TRUE, method = "REML", family = scat,
                na.action = na.exclude)
  # interaction term significant -->
  mco2te <- gam(co2corr ~ te(Time, DOY, bs = c("cc","tp")) + s(co2lag1),
                data = bdat, select = TRUE, method = "REML", family = scat,
                na.action = na.exclude)
  ## co2 dependence on previous value dwarfs all other terms... not sure I'm including everything
  ##    correctly? was i meant to do the autoregression separately?
  
  plot(acf(resid(mco2te), na.action = na.exclude)) # still autocorrelated.. could add AR2?
}
## try gamm approach following fromthebottomoftheheap?
## FIXME: should we pursue this?
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
testingnull <- gamm(co2corr ~te(Time, DOY, bs = c("cc","tp")), data = bdat, method="REML",
                    control = ctrl, verbosePQL = TRUE)
testing1 <- gamm(co2corr ~te(Time, DOY, bs = c("cc","tp")), data = bdat, method="REML",
                 correlation = corARMA(form = ~ 1, p=1),
                 control = ctrl, verbosePQL = TRUE)
##testing2 <- gamm(co2corr ~te(Time, DOY, bs = c("cc","tp")), data = bdat, 
       ##          correlation = corARMA(form = ~ 1, p=2),
      ##           control = ctrl, verbosePQL = TRUE)
## FIXME: "Error in `coef<-.corARMA`(`*tmp*`, value = value[parMap[, i]]) : 
##    Coefficient matrix not invertible""

res <- resid(testing1$lme, type = "normalized")
acf(res, lag.max = 200, main = "ACF - AR(2) errors")


## CO2 mechanistic model (2014): based on discussions with Helen, we should use rfu not mg/L
## ===========================================================================================
## no oxygen/pH here as it is taken to be one of the close covariables, without mechanistic info
##    per se (CO2, O2, pH)
mbio <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay), k=4) + 
                   s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
                   s(airtemp) + s(log(stability+1), k=3) + s(log(dailyrain+1), k=4) + s(conv) + 
                   s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
                 data = bdat, select = TRUE, method = "REML", family = tw,
                 na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
# had to put k for some vars to reduce overfitting
## FIXME: cond a bit strange!

if(runextras) {
  mbioti <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay), k=4) + 
                  s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
                  s(airtemp) + ti(log(stability+1)) + s(log(dailyrain+1), k=4) + ti(conv) +
                  s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
                data = bdat, select = TRUE, method = "REML", family = tw,
                na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
  
  mbiote <- gam(co2corr ~ s(chlrfu, by = factor(TimeofDay), k=4) + 
                  s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
                  s(airtemp) + ti(log(stability+1)) + s(log(dailyrain+1), k=4) + ti(conv) +
                  ti(log(stability +1), conv) +
                  s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
                data = bdat, select = TRUE, method = "REML", family = tw,
                na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
  anova(mbioti, mbiote, test="LRT")
  ## according to anova the interaction isn't important..
  ## FIXME -- gam's p-values have seemed a bit suspicious for a while... higher signif than expected
  ##    for many things... has something changed in computation? (e.g. windsp, with v low F)
}

## modeling co2 only based on oxygen and pH separated again since explains most everything
## tweedie was better than gamma based on AIC... not sure otherwise, both a bit off
if(runextras) {
  co2phmod <- gam(co2corr ~ s(ph1) + s(ODOrel1), method = "REML", family = tw(),
                  na.action = na.exclude, data=bdat, select=TRUE,
                  control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
  co2inter<- gam(co2corr ~ ti(ph1) + ti(ODOrel1) + ti(ph1, ODOrel1), method = "REML", family = tw(),
                 na.action = na.exclude, data=bdat,
                 control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
  ## FIXME -- why is the interaction significant when the relationship really just looks simple?
  anova(co2phmod, co2inter, test="LRT")
}

co2te <- gam(co2corr ~ te(ph1, ODOrel1), method = "REML", family = tw(),
             na.action = na.exclude, data=bdat, select=TRUE,
             control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))


## save models/output
## =========================================================================================

saveRDS(co2te, '../data/private/BPco2pho2mod.rds') # co2 te model with ph and o2
saveRDS(mbio, '../data/private/BPco2mechmod.rds') # mechanistic co2 model
saveRDS(testingnull, '../data/private/BPtimegammnull.rds') # time model no AR
saveRDS(testing1, '../data/private/BPtimegammAR1.rds') # time model AR1

## ==========================================================================================
## pH models (2014 and 2015) 
## 2015: no CO2 or O2 shallow data but other params ok
## ==========================================================================================


## Mechanistic models:
## ==========================================================================================
## Are the controls of pH the same in 2015?
## FIXME: note that some suspicious data in 2015, like unprecedented high cond, and high cdom..!??
##    has it been calibrated differently that year? Proceed to compare models with great caution
phbio <- gam(ph1 ~ s(log(chlrfu), by = factor(TimeofDay), k=4) + s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
               s(airtemp) + s(log(stability +1), k=3) + s(log(dailyrain+1), k=4) + s(conv) + 
               s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
             data = bdat, select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))
# had to put k for some vars to reduce overfitting

phbio5 <- gam(ph1 ~ s(log(chlrfu), by = factor(TimeofDay), k=4) + s(log(bga1rfu), by = factor(TimeofDay), k=4) + 
               s(airtemp) + s(log(stability +1), k=3) + s(log(dailyrain+1), k=4) + s(conv) + 
               s(log(turb), k=4) + s(cdom, k=3) + s(windsp) + s(cond1, k=4),
             data = bdat5, select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude, control = gam.control(nthreads = 3, trace = TRUE))

## Time models (2014, 2015):
## ==========================================================================================
## FIXME: using here the gamm approach, but if decided otherwise, need to mod
ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")

phgamm4AR1 <- gamm(ph1 ~te(Time, DOY, bs = c("cc","tp")), data = bdat, 
               correlation = corARMA(form = ~ 1, p=1), method="REML",
               control = ctrl, verbosePQL = TRUE)
phgamm5null <- gamm(ph1 ~te(Time, DOY, bs = c("cc","tp")), data = bdat5, method="REML",
                    control = ctrl, verbosePQL = TRUE)
phgamm5AR1 <- gamm(ph1 ~te(Time, DOY, bs = c("cc","tp")), data = bdat5, 
                correlation = corARMA(form = ~ 1, p=1), method="REML",
                control = ctrl, verbosePQL = TRUE) #for this model to run, I had to kill firefox, if I
#   could run it at all
#   "Error in recalc.corAR1(cSt, list(Xy = as.matrix(val))) : 
#   'Calloc' could not allocate memory (305795169 of 8 bytes)""

res <- resid(phgamm5AR1$lme, type = "normalized")
resnull <- resid(phgamm5null$lme, type = "normalized")

op <-par(mfrow=c(2,1))
acf(resnull, lag.max = 200, main = "ACF - default errors")
acf(res, lag.max = 200, main = "ACF - AR(1) errors")
par(op)

## trying bam to see if model fitting improves cause I can perhaps make more complex
## http://www.sfs.uni-tuebingen.de/~jvanrij/Tutorial/GAMM.html#setting-up-a-gamm-model
bamnull <- bam(ph1 ~te(Time, DOY, bs = c("cc","tp")), data = bdat5)

resnull <- resid(bamnull, type = "scaled.pearson")
valRho <- acf(resnull, plot=FALSE)$acf[2]

bamAR1 <- bam(ph1 ~te(Time, DOY, bs = c("cc","tp"), k=25), data = bdat5, rho=0.85, method='REML')
## FIXME: how to choose k intelligently? Can almost just choose to exactly follow the data
## NOTE! To get the AR1-accounted residuals here, need to use function 
acf_resid(bamAR1) #https://cran.r-project.org/web/packages/itsadug/vignettes/acf.html


## save best models
saveRDS(phbio, '../data/private/BPphmechmod2014.rds') # mecganistic pH model 2014
saveRDS(phbio5, "../data/private/BPphmechmod2015.rds") # mechanistic pH model 2015
saveRDS(phgamm5AR1, "../data/private/BPphtimegammAR12015.rds") # time pH model 2015 AR1
saveRDS(phgamm5null, "../data/private/BPphtimegammnull2015.rds") # time pH model 2015 no AR
saveRDS(phgamm4AR1, '../data/private/BPphtimegammAR12014.rds') # time pH model 2014 AR1
saveRDS(bamAR1, '../data/private/BPphtimebamAR12015.rds') # time pH model 2015 using BAM, AR1
