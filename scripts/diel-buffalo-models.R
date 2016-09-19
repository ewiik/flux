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

## gamming the time component of the data
m1 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(Week), data = bdat)

m2 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(Week) + ti(Time, Week, bs = c("cc","tp")), 
          data = bdat)
anova(m1, m2, test = "LRT") # m2 better

m3 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")), 
          data = bdat)
m4 <- gam(ph1 ~ te(Time, DOY, bs = c("cc","tp")), data=bdat)

mco2 <- gam(co2corr ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")), 
            data = bdat)

plot(m3, scheme=2) 
plot(acf(resid(m3)))

pdf("diel-bpgam.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(m2, pages = 4, scheme = 2)
par(op)
dev.off()

saveRDS(m4, "../data/private/bp-diel-timemod2014.rds")


## gam on other params
## ===========================================================================================
bdat <- transform(bdat, daylength = sunset-sunrise)
bdat$daylength <- as.numeric(bdat$daylength)
bpmod <- gam(ph1 ~
               s(chl, k=3) + #was overfitting
               s(airtemp, k=3) +
               s(daylength, k=4) +  #was overfitting
               #s(ODOrel1, k=4),
               s(ODOabs1), 
             data = bdat,
             select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE,
                                   newton = list(maxHalf = 60)))
## FIXME: add time of day
## FIXME: tested scat but though histogram looked better nothing else did...??
## tested windsp, meanwindsp, and lag1.. nothing amazing
## odo rather close to strat so dropped strat here
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
