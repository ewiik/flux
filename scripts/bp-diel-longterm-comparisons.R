## compare diel and longterm for bp

## load packages
library('ggplot2')
library('mgcv')
library('viridis')

## read in data
blong <- readRDS('../data/private/bp-longterm-mod.rds')
bdiel <-  readRDS('../data/private/bpbuoy2014-mod.rds')

## read in model
source("../functions/gasExchangeUser.R")

## subset data
bdiel <- subset(bdiel, Hour >= 10 & Hour <=15)
bdiel$Date <- as.Date(bdiel$datetime)

bdat <- subset(blong, year == 2014, select = c('ph', 'CO2fluxmmolm2d', 'pCO2uatm', 'week','year'))
isna <- which(is.na(bdat$pCO2uatm))
bdat <- bdat[-isna,]

## try to compare data by week...
bdatall <- merge(bdat, bdiel, by.x = 'week', by.y = 'Week')
with(bdatall, plot(pCO2uatm ~ co2corr)) #nooooo

## create all Mondays for time series ffs
allMondays <- as.Date(seq(c(ISOdate(1984,12,31)), by="7 DSTdays", length.out=1700))
allMondays <- data.frame(Date = allMondays)
allMondays$Week <- format(allMondays$Date, '%V')
allMondays$Year <- format(allMondays$Date, '%Y')

## merge
bdatall <- merge(bdat, allMondays, by.x = c('week','year'), by.y = c('Week','Year'))

## take all days, and only days that are the same
with(bdiel,plot(co2corr ~ Date, col='red', ylab='pCO2', ylim=c(0,2000)))
points(bdatall$pCO2uatm ~ bdatall$Date, pch=21, bg='black')
legend('topright', legend = 'red=measured \n black=calculated', border = FALSE)

bmatch <- which(bdiel$Date %in% bdatall$Date)
bmatched <- bdiel[bmatch,]
with(bmatched,plot(co2corr ~ Date))
points(bdatall$pCO2uatm ~ bdatall$Date, pch=3)
legend('bottomleft', legend = 'o = meas \n+ = calc', border = FALSE)

## long term variation in relationship between dic and cond
ggplot(blong, aes(x = cond, y = totaldic, group=year, colour=year)) +
  geom_point() +
  scale_color_viridis()
  #facet_wrap('year')

blong$year <- as.factor(blong$year)
dicmod <- gam(totaldic ~ s(cond) + s(cond, by = year, m = 1, k=3) + s(year, bs = "re"), data=blong,
              method = "REML", family = gaussian,
              na.action = na.exclude,
              control = gam.control(nthreads = 3, trace = TRUE))
dicmodnull <-gam(totaldic ~ s(cond) + s(year, bs = "re"), data=blong,
                 method = "REML", family = gaussian,
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE))
anova(dicmodnull, dicmod, test = 'LRT') # --> more complex model warranted

newdat <- data.frame(year = rep(2014), cond = bdiel$cond1)
preddic <- predict(dicmod, newdata=newdat, type = 'response')
newdat$dic <- preddic
with(newdat, plot(dic ~ cond))

bdiel <- transform(bdiel, kerridic = 26.57 + 0.018 * cond1)  # (mg/L)
bdiel$kerridicumol <- bdiel$kerridic / 0.012
plot(newdat$dic ~ bdiel$kerridic)
abline(0,1)
bdiel$gamdicumol <- newdat$dic / 0.012
pco2kerri <- with(bdiel, gasExchangeUser(temp=temp1, cond = cond1, dic=kerridicumol, salt=NULL,ph=ph1, 
                                         wind=windsp, alknotdic = FALSE, 
                                         kpa=pressureKPA, pco2atm = maunaloa))
pco2gam <- with(bdiel, gasExchangeUser(temp=temp1, cond = cond1, dic=gamdicumol, salt=NULL,ph=ph1, 
                                         wind=windsp, alknotdic = FALSE, 
                                         kpa=pressureKPA, pco2atm = maunaloa))
bdiel$pco2gam <- pco2gam$pco2
bdiel$pco2kerri <- pco2kerri$pco2

with(bdiel, plot(pco2gam ~ pco2kerri))
abline(0,1)

ggplot(bdiel, aes(x=pco2kerri, y=pco2gam)) + geom_point() + geom_abline(slope=1, intercept=0)
ggplot(bdiel, aes(x=pco2kerri, y=ph1)) + geom_point() 
ggplot(bdiel, aes(x=pco2gam, y=ph1)) + geom_point() 

