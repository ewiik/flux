## statistics on bp data

## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  source("diel-buffalo.R")}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')
bdat <- bdat[order(bdat$datetime),]

bdat5 <- readRDS('../data/private/bpbuoy2015-mod.rds')
bdat5 <- bdat5[order(bdat5$datetime),]

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
lmnight <- with(bdat[bdat$isDay == FALSE,], lm(co2corr~pco2))
lmday <- with(bdat[bdat$isDay == TRUE,], lm(co2corr~pco2))
lmint <- with(bdat, lm(co2corr~pco2))
with(bdat, cor(x=co2corr, y=pco2, use='complete.obs', method='pearson'))

# just with the 10am-3pm data
ten3 <- subset(bdat, Hour >= 10 & Hour <=15)

## could do this to make converge but also did manually, required four iterations
reps <- 8 #takes a while for it to converge
for (i in 1:reps) {
   start <- if (i == 1 ) {coef(glm(co2corr ~ pco2, data = ten3, family = gaussian))}
     else { coef(glmmod) }
   glmmod <- glm(co2corr ~ pco2, data = ten3, family = Gamma(link = "identity"), start = start,
                 control = glm.control(maxit=100))
 }
## FIXME: this is not an adequate model based on QQ plot so need to rethink at some point

# residual plot etc
plot(resid(lmint) ~ lmint$fitted.values, ylab='Residuals (pCO2 uatm)', 
     xlab='Fitted values (pCO2 uatm)')

## gamming the time component of the data for pH
## =================================================================================================
## create lag terms
bdat <- transform(bdat, timelag1 = myAR(bdat$Time, 1))

## create models; testing AR inclusion and comparing k with edf to set k (gam.check)
mnull <- gam(ph1 ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=30), data=bdat) 
m1 <- gam(ph1 ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=30) + ti(Time, DOY, bs = c("cc","tp")), 
          data = bdat)
m2 <- gam(ph1 ~ te(Time, DOY, bs = c("cc","tp")), data=bdat)

plot(acf(resid(m1)))

anova(mnull, m1, test = "LRT") # m1 better

## apply same to co2
mco2null <- gam(co2corr ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=30), data=bdat, 
                select = TRUE, 
                method = "REML", family = scat, na.action = na.exclude)
# + s(timelag1, bs='cc'), 
mco2ti <- gam(co2corr ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=30) + s(timelag1, bs='cc') + 
              ti(Time, DOY, bs = c("cc","tp"), k=15), 
            data = bdat, select = TRUE, method = "REML", family = scat,
            na.action = na.exclude) ## crashes with higher k for interaction!
anova(mco2null, mco2ti, test = "LRT") # mco2 better

mco2te <- gam(co2corr ~ s(timelag1, bc='cc') + 
              te(Time, DOY, bs = c("cc","tp")), # worried of setting k higher
            data = bdat, select = TRUE, method = "REML", family = scat,
            na.action = na.exclude) 

plot(acf(resid(mco2null))) ## FIXME: is < 0.2 ok?
plot(mco2, scheme=2) ## note that splines are almost flat for DOY and Time for all models
##    with AR1 term

mco2full <- gam(co2corr ~ s(Time, bs = "cc", k=20) + s(DOY, k=30) + s(Week) + 
                  s(co2lag1), data=bdat, 
                select = TRUE, 
                method = "REML", family = scat, na.action = na.exclude)

pdf("../docs/private/diel-bpgam.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(m2, pages = 4, scheme = 2)
par(op)
dev.off()

saveRDS(m4, "../data/private/bp-diel-timemod2014.rds")

testing <- gam(co2corr ~ ti(Time, bs = "cc", k=20) + ti(DOY, k=30) + s(timelag1, bs='cc', k=20), 
               data=bdat, 
                select = TRUE, 
                method = "REML", family = scat, na.action = na.exclude)

ctrl <- list(niterEM = 0, msVerbose = TRUE, optimMethod="L-BFGS-B")
testing2 <- gamm(co2corr ~ s(Time, bs = "cc", k = 20) + s(DOY, k = 30),
                    data = bdat, correlation = corCAR1(form = ~ datetime), #Time|DOY
                      control = ctrl)
## gamm crashes if I specify 0.8 as the correlation!
res <- resid(testing2$lme, type = "normalized")
acf(res, lag.max = 36, main = "ACF - AR(2) errors")

## gam on other params
## ===========================================================================================
bpmod <- gam(co2corr ~
               s(chl, k=3) + #was overfitting
               s(bga1cell) +
               s(airtemp, k=3) +
               s(schmidt) +  
               s(ODOrel1, k=4) +
               s(ph1)  +
               s(ODOabs1), 
             data = bdat,
             select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE,
                                   newton = list(maxHalf = 60)))
## FIXME: add time of day
## FIXME: tested scat but though histogram looked better nothing else did...??
plot(bpmod, pages=1)
saveRDS(bpmod, '../data/private/bpmod.rds')

bpco2phmod <- gam(co2corr ~ s(ph1), method = "REML", family = gaussian,
                  na.action = na.exclude, data=bdat2014fullf,
                  control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
saveRDS(bpco2phmod, '../data/private/bpco2phmod.rds')

## can temporary stratificatino be predicted with wind and what is the best lag?
bdat <- transform(bdat, strat=temp2-temp1)
stratmod <- gam(strat ~ s(windsp) + s(winddir), data = bdat,
                select = TRUE, method = "REML", family = gaussian(),
                na.action = na.exclude)
windsplit <- with(bdat, split(bdat, list(Month, Day)))
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(Month = df['Month'][1,], Day = df['Day'][1,]), t(means))
}
windmeans <- do.call(rbind, lapply(windsplit, mycolMeans, cols=c("windsp", "winddir")))
windmeans <- windmeans[order(windmeans$Month, windmeans$Day),]
windmeans <- na.omit(windmeans)
names(windmeans)[grep("wind", names(windmeans))] <- c("meanwindsp", "meanwinddir")
for(i in 1:nrow(windmeans)) {
  if(i == 1) { windmeans$lag1[i] <- NA} else {
  windmeans$lag1[i] <- windmeans$meanwindsp[i-1]
  }}
for(i in 1:nrow(windmeans)) {
  if(i <= 2) { windmeans$lag2[i] <- NA} else {
    windmeans$lag2[i] <- windmeans$meanwindsp[i-2]
  }}

bdat <- merge(bdat, windmeans)
stratmod <- gam(strat ~ s(meanwindsp, k=3) + s(lag1, k=3) + s(temp1, k=3) + s(meanwinddir), data = bdat,
                select = TRUE, method = "REML", family = gaussian(),
                na.action = na.exclude)

## a model for oxygen?
oxmod <- gam(ODOrel1 ~
               s(chl, k=3) + #was overfitting
               s(airtemp) +
               s(daylength, k=4) +  #was overfitting
               s(strat, k=3),
             data = bdat,
             select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE,
                                   newton = list(maxHalf = 60)))
plot(oxmod, pages = 1)

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
