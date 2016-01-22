## having done the CO2 flux ~ ... regression for high-reso routines data, this script
##    carries on and investigates what may be driving (covarying with) pH
## We know that the max that pH can achieve is related to the start-of-year pH, so
##    therefore we will also look at Year as a random effect!
## We know that N and P tend to be correlated, so let's allow one or both to be dropped
##    and start with allowing shrinkage

## get necessary packages
library("ggplot2")
library("reshape2")
library("mgcv")

## load the data we need
if (!file.exists("data/private/regvars.rds")) {
  source("scripts/regression_routines.R")
}
regvars <- readRDS("data/private/regvars.rds")

## is pH heavy-tailed?
plot(density(regvars$pH_surface, na.rm = TRUE))
qqnorm(regvars$pH_surface, na.rm = TRUE)
## yes, choose family = scat() in gam

## create model; pH a function of instantaneous productivity and previous year sum index?
## FIXME: thought I should transform Year into factor since it's a random effect, but
##    for some reason actually model didn't work with Year as not a factor,...
##    (as s(Year) or s(Year, bs = "re")) which is weird
regvarf <- regvars
regvarf <- transform(regvarf, Year = as.factor(Year))
phmod <- gam(pH_surface ~ s(Lake, bs = "re") + s(Year, bs = "re") + s(Chl_a_ug_L) + s(GPP_h) + 
               s(TDN_ug_L) + s(DOC_mg_L) + s(Oxygen_ppm) + te(PDO, SOI), data = regvarf,
             select = TRUE, method = "REML", family = scat()) 
summary(phmod)
plot(phmod, pages = 1, pers = TRUE)

## A miniscript that will run the plots for each individual variable's effect keeping
##    all other variables constant (at their mean)
## Note that ..$model is the data that was used in the modeling process
varmeans <- data.frame(t(colMeans(phmod$model[,4:10])))
varlist <- names(varmeans)
varlist
varying <- "DOC_mg_L"
varindexv <- which(names(varmeans) == varying)
varindexm <- which(names(phmod$model) == varying)

testdata <- data.frame(cbind(varmeans[,-varindexv], phmod$model[,varindexm]))
names(testdata)[length(testdata)] <- varying

testdata$pH_surface <- phmod$model$pH_surface
testdata$Lake <- phmod$model$Lake
testdata$Year <- phmod$model$Year
fits  <-  predict(phmod, newdata=testdata, type='response', se = TRUE)
predicts  <-  data.frame(testdata, fits)
names(predicts)[names(predicts) == varying]
fit <- "fit"
ggplot(aes_string(x=varying,y=fit), data=predicts) + # aes_string means it takes the colname
  #   indicated by the string that varying refers to!!
  geom_smooth(aes(ymin = fit - 1.96*se.fit, ymax=fit + 1.96*se.fit),
              fill='gray80', size=1,stat='identity') +
  xlab(varying)
## FIXME: discuss plots with Gavin, especially the zigzag ones

## create crude correlate of previous year's autumn productivity to see if signif
##    using chl a
## grab last two possible chl a data points
lastprod <- subset(regvars, Month == 8 | Month == 9, select = c(Lake, Year, Date, Chl_a_ug_L))
lastspl <- with(lastprod, split(lastprod, list(Lake, Year)))

## create wrappers that do colmeans for the rows we want and spits out df we can merge
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(Lake = df['Lake'][1,], Year = df['Year'][1,]), t(means))
}

grablast <- function(df) {
  grab <- c(nrow(df), nrow(df) -1)
  df <- df[grab,]
}

## apply wrappers
lastspl <- lapply(lastspl, grablast)
chlmeans <- lapply(lastspl, mycolMeans, cols = c("Chl_a_ug_L")) 
chlmeans <- do.call(rbind, chlmeans)
rownames(chlmeans) <- NULL

## tidy chlmeans for merging: make mean apply to following year and use that to merge
chlmeans <- transform(chlmeans, PrevYear = Year + 1)
chlmeans <- chlmeans[,-which(names(chlmeans) == 'Date' | names(chlmeans) == 'Year')]
regvarl <- regvars

## merge by appropriate year
regvarl <- merge(regvarl, chlmeans, by.x = c('Lake', 'Year'), by.y = c('Lake', 'PrevYear'))
regvarl <- transform(regvarl, Year = as.factor(Year))

## create new regression to see if any sense in including it
phmod2 <- gam(pH_surface ~ s(Lake, bs = "re") + s(Year, bs = "re") + s(Chl_a_ug_L.x) + 
                s(Chl_a_ug_L.y) +s(GPP_h) + s(TDN_ug_L) + s(DOC_mg_L) + s(Oxygen_ppm) + 
                te(PDO, SOI), data = regvarl,
             select = TRUE, method = "ML", family = scat()) 
summary(phmod2)
## nope, no sense

