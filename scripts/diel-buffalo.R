## load in Helen's 2015 buoy data -- from file sent 28.01.2016
## And 2014 data from email sent Monday - April 25, 2016 
## two sonde depths, n.1 is the shallower one based on PAR..
## Not sure what cdom and mv are measures of.. and why so many temps?
## FIXME: really need to learn what's what with the column name matching (2015)
## FIXME: removed outliers but also 0s that are not real 0s.. need to remove (2015)

## load necessary packages
library("ggplot2")
library('StreamMetabolism')
library('mgcv')
library('viridis')
library('reshape2')

## load necessary functions
source('../functions/gasExchangeUser.R')
source('../data/private/salcalc.R')

## load in supplementary data
regvars <- readRDS('../data/private/regvars.rds')
params <- readRDS('../data/private/params-flux.rds')
buffreg <- subset(regvars, select = c('Date', 'Lake', 'Month', 'DOY', 'pH_surface',
                                      'lakepCO2','co2Flux'), Lake == 'B' & Year > 2013)
## mauna loa for 2014 months, such few data points I copy pasted manually from
##    ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt
ml <- read.csv('../data/maunaloa2014summer.csv')
names(ml)[which(names(ml)=='pCO2')] <- 'maunaloa'
ml <- ml[,-which(names(ml) == 'pCO2interp')]

## ======================================================================================
## 2014
## ======================================================================================

## read in data; initially got only part of the supporting info, but the co2 that Helen
##    already corrected. Then asked for pH etc and she sent that but it doesn't have
##    the corrected values
bdat2014 <- read.csv('../data/private/BuffPd_DataBlobs_Helen.csv', skip = 4, 
                     col.names = c('datetime', "winddir", "windsp", "airtemp", 
                                   "relhum", "pressure", "dailyrain", "par1", "par2",
                                   "par3", "temp1", "ODOrel1", "ODOabs1", "temp2",
                                   "temp3", "temp4", "temp5", "co2corr"))
bdat2014supp <- read.csv("../data/private/BPBuoy2014.csv", skip = 4, 
                         col.names = c('datetime',"batvolt", "winddir", "windsp", "airtemp", 
                                       "relhum", "pressure", "dailyrain", "par1", "par2",
                                       "par3", "cdom", "mvolts", "co2-1", "co2-2", "temp1",
                                       "cond1", "ph1", "ph1mV", "turb", "chl", "chlrfu", "bga1cell",
                                       "bga1rfu", "ODOrel1", "ODOabs1", "temp2", "cond2", "ph2", 
                                       "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3", 
                                       "temp4", "temp5", "temp6"))

## make time understandable
bdat2014 <- transform(bdat2014, datetime = as.POSIXct(as.character(datetime), 
                                              format = "%m/%d/%Y %H:%M"))
bdat2014 <- transform(bdat2014, Hour = as.numeric(format(datetime, format = "%H")))
bdat2014 <- transform(bdat2014, Month = as.numeric(format(datetime, format = "%m")))
bdat2014 <- transform(bdat2014, Day = as.numeric(format(datetime, format = "%d")))
bdat2014 <- transform(bdat2014, DOY = as.numeric(format(datetime, format = "%j")))
bdat2014 <- transform(bdat2014, Time = strftime(datetime, format = "%H:%M", tz="GMT-1"))
bdat2014 <- transform(bdat2014, isDay = ifelse(Hour < 21 & Hour > 7, TRUE, FALSE))

bdat2014supp <- transform(bdat2014supp, datetime = as.POSIXct(as.character(datetime), 
                                                      format = "%m/%d/%Y %H:%M"))
bdat2014supp <- transform(bdat2014supp, Hour = as.numeric(format(datetime, format = "%H")))
bdat2014supp <- transform(bdat2014supp, Month = as.numeric(format(datetime, format = "%m")))
bdat2014supp <- transform(bdat2014supp, Day = as.numeric(format(datetime, format = "%d")))
bdat2014supp <- transform(bdat2014supp, DOY = as.numeric(format(datetime, format = "%j")))
bdat2014supp <- transform(bdat2014supp, Time = as.numeric(format(datetime, "%H")) +
  as.numeric(format(datetime, "%M"))/60)
bdat2014supp <- transform(bdat2014supp, Week = as.numeric(format(datetime, "%U")))
bdat2014supp <- transform(bdat2014supp, isDay = ifelse(Hour < 21 & Hour > 7, TRUE, FALSE))

## in 2014 pressure is in HPA but need KPA for CO2 calcs
bdat2014supp$pressureKPA <- bdat2014supp$pressure/10

## cleaning dates in 2014 (sensor not cleaned in 2015)
dates <- c('12/08/2014', '15/07/2014', '19/07/2014', '26/06/2014') # last may not be a cleaning
#   date
dates <- as.POSIXct(dates, format = "%d/%m/%Y")
datesDOY <- format(dates, "%j")

## homogenise data
bdat2014supp$pressureKPAwater <- bdat2014supp$pressureKPA + .08
bdat2014supp <- merge(bdat2014supp, ml)

## calculate salinity
bdat2014supp <- transform(bdat2014supp, salcalc = salcalc(temp=temp1, cond=cond1, 
                                                         dbar=pressureKPAwater/10))
bdat2014suppflux <- with(bdat2014supp, gasExchangeUser(temp = temp1, cond = cond1, ph=ph1, 
                                                       wind=windsp, alknotdic = FALSE, 
                                                       salt = salcalc, kpa=pressureKPA,
                                                       pco2atm = maunaloa, diccalc = TRUE))
bdat2014suppall <- cbind(bdat2014supp, bdat2014suppflux)
bdat2014corr <- subset(bdat2014, select=c('datetime', 'co2corr'))
bdat2014full <- merge(bdat2014suppall, bdat2014corr, all.x = TRUE)

## let's look at when sun up and down
sundat <- sunrise.set(lat = 50.648016, long=-105.5072930, date = '2014/06/01', 
                      timezone = 'UTC+6', num.days=90)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

## merge sundat with BP data
bdat2014full <- merge(bdat2014full, sundat)

## change isday to depend on sunup and down rather than rigid times
bdat2014full$isDay <- ifelse(bdat2014full$Time < bdat2014full$DownTime & 
                               bdat2014full$Time > bdat2014full$UpTime, 
                             TRUE, FALSE)
bdat2014full$TimeofDay <- ifelse(bdat2014full$Time < bdat2014full$DownTime & 
                               bdat2014full$Time > bdat2014full$UpTime, 
                             'Day', 'Night')


## save object for other purposes
saveRDS(bdat2014full, '../data/private/bpbuoy2014-mod.rds')

## ==========================================================================================
## 2015
## ==========================================================================================
bdat <- read.csv("../data/private/BPBuoyData2015raw.csv", skip = 4, 
                 col.names = c("datetime", "batvolt", "winddir", "windsp", "airtemp", 
                               "relhum", "pressure", "dailyrain", "par1", "par2",
                               "par3", "cdom", "mvolts", "co2-1", "co2-2", "temp1",
                               "cond1", "ph1", "ph1mV", "turb", "chl", "chlrfu", "bga1cell",
                               "bga1rfu", "ODOrel1", "ODOabs1", "temp2", "cond2", "ph2", 
                               "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3", 
                               "temp4", "temp5", "temp6", "temp7"))

press <- read.csv("../data/bp-pressure2015.csv")
names(press)[which(names(press) == 'Pressure')] <- 'AirPressure'
press <- press[,-which(names(press) == 'Year')]

bdat <- transform(bdat, datetime = as.POSIXct(as.character(datetime), 
                                              format = "%m/%d/%Y %I:%M:%S %p"))
bdat <- transform(bdat, Hour = as.numeric(format(datetime, format = "%H")))
bdat <- transform(bdat, Month = as.numeric(format(datetime, format = "%m")))
bdat <- transform(bdat, Day = as.numeric(format(datetime, format = "%d")))
bdat <- transform(bdat, DOY = as.numeric(format(datetime, format = "%j")))

bdat <- transform(bdat, isDay = ifelse(Hour < 21 & Hour > 7, TRUE, FALSE))

## get pressure in
bdat <- merge(bdat, press, by=c('Month', 'Day'))

## check the CO2 columns
with(bdat, plot(co2.1 ~ datetime))

## not sure why all june and july values are completely different for both co2 profiles..!?     
with(bdat, plot(density(co2.1), main='Density of 2015 uncorrected all-outliers-in pCO2'))
with(bdat, plot(density(co2.2)))

## remove the crazy negative and positive values and bound the ppm
##    but first make note of which rows these are
outliers <- which(abs(bdat$co2.1) > 1200 | abs(bdat$co2.2) > 1200)
bdat<- bdat[-outliers,]
## still bimodal but don't wanna remove too much in case it's real

## initial pressure in hg inches or kPa (met), need to get hPa for corrections
bdat$pressure[grep('-Invalid-', bdat$pressure)] <- NA
bdat$pressure <- as.character(bdat$pressure)
bdat$pressure <- as.numeric(bdat$pressure)
bdat$pressureHPA <- bdat$pressure/0.0295299830714

bdat$AirPressureHPA <- bdat$AirPressure * 10

## deal with similar issues with wind (is in m/s)
bdat$windsp[grep('-Invalid-', bdat$windsp)] <- NA
bdat$windsp <- as.character(bdat$windsp)
bdat$windsp <- as.numeric(bdat$windsp)

## plot some diagnostics and look if pH and oxygen determining CO2
##    assuming I know which odo and ph measure relates to which co2 probe...
with(bdat, plot(co2.1 ~ Hour, col = ifelse(isDay, "red", "black")))
with(bdat, plot(co2.2 ~ Hour, col = ifelse(isDay, "red", "black")))

with(bdat, plot(co2.1 ~ ph1, col = ifelse(Month < 9, "red", "black")))
with(bdat, plot(co2.1 ~ ODOrel1, col = ifelse(isDay, "red", "black")))
## something wrong with ODO measurements!! loads 0s. Let's remove 0s
## Existing data mostly contiguous in time
with(bdat[-which(bdat$ODOrel1 == 0), ], plot(co2.1 ~ ODOrel1, col = ifelse(Day, "red", "black")))


with(bdat, plot(co2.2 ~ ph2, col = ifelse(Month < 9, "red", "black")))
with(bdat, plot(co2.2 ~ ODOrel2, col = ifelse(Month < 9, "red", "black")))
## here, something fishy with the month of September.. CO2 mostly 0

with(bdat[-which(bdat$ODOrel1 == 0), ], 
     plot(ODOrel1 ~ ph1, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
with(bdat, plot(ODOrel2 ~ ph2, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
## night and day different slopes in sonde1, doesn't happen in sonde2 even when
##    following the ODOrel1 subset

with(bdat[-which(bdat$ODOrel1 == 0), ], 
     plot(ODOrel1 ~ chl, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
with(bdat, plot(ph1 ~ chl, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))


## correction as per ecohydrology paper
#The post-measurement correction of sensor output as pCO2 related to changes in temperature 
#and pressure were determined empirically for the Vaisala GMT IRGAs (Vaisala Oyj, 2008). 
#Pre-corrected sensor output needs to be reduced by 0.15% of the measured reading per hPa 
#increase in pressure relative to calibration pressure (typically 1013 hPa). Pressure readings
#below the calibration pressure require increasing the sensor output by 0.15% of the measured
#reading per hPa. Pre-corrected sensor output also needs to be increased by 0.3% of the measured 
#reading per ?C of increased temperature relative to calibration temperature (typically 25?C). 
#Water temperature above the calibration temperature requires decreasing the sensor output by 
#0.3% of the measured reading per ?C. An additional correction is required when sensors are
#deployed in an aquatic environment. Water depth affects the pressure exerted on the sensor, 
#and this water depth correction is added to the atmospheric pressure correction. For example, 
#10 cm water depth above the sensor location corresponds to an increase in pressure of 9.81 hPa. 
#If sensor output is 1000 ppm pCO2, the output needs to be reduced by 14.72 ppm because of 
#the additional pressure exerted by the increased water level. If the water depth varies 
#relative to the sensor location, a depth correction must be determined for each measurement 
#time interval (e.g. by recording water depth using a pressure transducer or other water level sensor)
corrpress <- bdat$co2.1*((1013 - bdat$AirPressureHPA) * 0.0015)
corrtemp <- bdat$co2.1*((25 - bdat$temp1) * 0.003)
corrairpress <- -bdat$co2.1*((1014.14 - 1013.25) * 0.0015)
# 89cm pressure calculated at http://www.calctool.org/CALC/other/games/depth_press
#   is where the CO2 sensor is.

bdat$co2.1corr <- bdat$co2.1 + corrpress + corrtemp + corrairpress

## plot pH with co2 
bdat2 <- bdat
bdat2$Month <- as.factor(bdat2$Month)
ggplot(data=bdat2, aes(y = co2.1corr, x = ph1, col = Month)) +
  geom_point() +
  ylab('CO2 (ppm) corrected') +
  xlab('pH') +
  geom_text(label='2015', x=8.9, y=1250, col='black')

## calculate CO2 for 2015
bdat$dic <- rep(NA)
bdatflux <- with(bdat, gasExchangeExtra(temp = temp1, cond = cond1, ph=ph1, wind=windsp, alknotdic = FALSE,
                                       kerri = TRUE, dic = dic, kpa=AirPressure))
bdatall <- cbind(bdat, bdatflux)

ggplot(data=bdatall, aes(y = co2.1corr, x = pco2, col = factor(Month))) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  geom_abline(intercept=0,slope=1) +
  ylab("pCO2 (measured)") +
  xlab("pCO2 (calculated)") +
  geom_text(label='2015', x=200, y=1250, col='black')

lm2015 <- with(bdatall, lm(co2.1corr~pco2 -1))
with(bdatall, cor(co2.1corr, pco2,use='complete.obs', method='pearson'))

## plot 2015 data
ggplot(data=bdat, aes(y = co2.1, x = DOY, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  ylab("CO2 (ppm)") +
  xlab("Day of year 2015")

