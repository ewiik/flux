## load in Helen's 2015 buoy data -- from file sent 28.01.2016
##    UPDATE: now checked and commented to be consistent with file sent by Julie 15thSep2016!
## And 2014 data from email sent Monday - April 25, 2016 
##    UPDATE: now checked and commented to be consistent with file sent by Julie 15thSep2016!
##    However note that co2.1 goes to Sep and co2corr not so this would be nice to amend.

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
## Note: "pHDeep showing suspicious values around the time of a cleaning by the Regina 
##    group. pH is < 8 (for the only time during the season) just before the cleaning 
##    and then increases to > 9 just afterwards." (notes file received from Julie)
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
testing2014 <- read.csv("../data/private/BPoundBuoy_2014_QC_Processed_Data.csv") # most up to date 
#   version received from Julie on Sep 15th 2016

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

testing2014 <- transform(testing2014, DateTime = as.POSIXct(as.character(DateTime), 
                                                    format = "%Y-%m-%d %H:%M"))
testing2014 <- transform(testing2014, Hour = as.numeric(format(DateTime, format = "%H")))
testing2014 <- transform(testing2014, Month = as.numeric(format(DateTime, format = "%m")))
testing2014 <- transform(testing2014, Day = as.numeric(format(DateTime, format = "%d")))
testing2014 <- transform(testing2014, DOY = as.numeric(format(DateTime, format = "%j")))
testing2014 <- transform(testing2014, Time = as.numeric(format(DateTime, "%H")) +
                       as.numeric(format(DateTime, "%M"))/60)

missing <- which(!testing2014$DateTime %in% bdat2014supp$datetime )
testing2014[missing,1]## these are the confusing row shift start columns where co2corr et al are
##    mismatched between supp and corr file, and 3rd-4th Sep
missingalt <- which(!testing2014$DateTime %in% bdat2014$datetime )
testing2014[missingalt,1]
## hmm many more missing here; those above but the latter starting from 22.08.

missing2 <- which(!bdat2014supp$datetime %in% testing2014$DateTime)
missing2alt <- which(!bdat2014$datetime %in% testing2014$DateTime)
## none here. Basically the newest master file has MORE dates than the ones I was previously
##    given.... notes say that some sondes quit on 22nd Aug but data continued to 4th Sep.
## Can live without the 3rd and 4th of Sep, and grabbed co2corr no others from corr, so will 
##    proceed with checking the start of the series
head(testing2014$WindDir, n= 10) # 19:30 is 214, so supp is correct and corr has been corrupted
#[1] 269 198 242 220 223 214 200 216 202  79
head(bdat2014$winddir, n= 10)
#[1] 198 242 220 223 214 200 216 202  79  68
head(bdat2014supp$winddir, n= 10)
#[1] 214 200 216 202  79  68  74  80  81  76

## I need to subtract 1h from all vars in corr!! ugh   
bdat2014$datetime <- bdat2014$datetime - 60*60

## now remove trailing septembers from testing, and the up to 19:30 at the start
testing2014 <- testing2014[-missing,]

## testing has some -10000 where supp has NA..
odds <- which(testing2014$CO2shallow <= 0)
testing2014$CO2shallow[odds] <- NA

## need to remove the start ones from bdat2014 after the hourshift
missing2alt2 <- which(!bdat2014$datetime %in% testing2014$DateTime)
bdat2014 <- bdat2014[-missing2alt2,]

## ok now we have complete match of dates... let's check other vars!
plot(bdat2014supp$co2.1 ~ testing2014$CO2shallow, type = "l")
plot(bdat2014supp$ph1 ~ testing2014$pHShallow, type = "l")
plot(bdat2014supp$windsp ~ testing2014$WindSp, type = "l")

## and let's check that the corr now is in agreement
plot(bdat2014$winddir ~ 
       testing2014$WindDir[testing2014$DateTime <= as.POSIXct("22.08.2014 07:45", 
                                                              format="%d.%m.%Y %H:%M")], 
     type = "l")

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

## yep that is good;now for flags that I might cause me to take more rows out or make NA
## 1. maintenance/calibration
bdat2014full[bdat2014full$datetime >= as.POSIXct("26.06.2014 10:45", format="%d.%m.%Y %H:%M") &
       bdat2014full$datetime <= as.POSIXct("26.06.2014 13:00", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat2014full[bdat2014full$datetime >= as.POSIXct("08.07.2014 09:45", format="%d.%m.%Y %H:%M") &
       bdat2014full$datetime <= as.POSIXct("08.07.2014 17:45", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat2014full[bdat2014full$datetime >= as.POSIXct("29.07.2014 10:45", format="%d.%m.%Y %H:%M") &
       bdat2014full$datetime <= as.POSIXct("29.07.2014 11:30", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat2014full[bdat2014full$datetime >= as.POSIXct("31.07.2014 09:45", format="%d.%m.%Y %H:%M") &
       bdat2014full$datetime <= as.POSIXct("31.07.2014 17:15", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat2014full[bdat2014full$datetime >= as.POSIXct("12.08.2014 10:15", format="%d.%m.%Y %H:%M") &
       bdat2014full$datetime <= as.POSIXct("12.08.2014 11:00", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat2014full[bdat2014full$datetime >= as.POSIXct("14.08.2014 11:00", format="%d.%m.%Y %H:%M") &
       bdat2014full$datetime <= as.POSIXct("14.08.2014 21:15", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat2014full[bdat2014full$datetime >= as.POSIXct("26.08.2014 08:15", format="%d.%m.%Y %H:%M") &
               bdat2014full$datetime <= as.POSIXct("26.08.2014 15:30", format="%d.%m.%Y %H:%M"),-1] <- NA
## 2. General failure
bdat2014full[bdat2014full$datetime >= as.POSIXct("22.08.2014 08:00", format="%d.%m.%Y %H:%M") &
               bdat2014full$datetime <= as.POSIXct("22.08.2014 11:00", format="%d.%m.%Y %H:%M"),-1] <- NA

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

## remove outliers
bdat2014full$cdom[bdat2014full$cdom < 0]  <- NA
bdat2014full$cdom[bdat2014full$cdom < 40 & bdat2014full$Month >= 8] <-NA
bdat2014full$cdom[bdat2014full$cdom < 43.5 & bdat2014full$Month >= 8 & bdat2014full$Day > 3] <-NA

bdat2014full$turb[log10(bdat2014full$turb) > 2.9] <- NA

bdat2014full$par2[bdat2014full$par2 < 0] <- NA # deal with -100 000 values
bdat2014full$par3[bdat2014full$par3 < 0] <- NA

## save object for other purposes
saveRDS(bdat2014full, '../data/private/bpbuoy2014-mod.rds')

## ==========================================================================================
## 2015: only using surface since deep just kept malfunctioning all year
## ==========================================================================================
bdat <- read.csv("../data/private/BPBuoyData2015raw.csv", skip = 3, 
                 col.names = c("datetime", "batvolt", "winddir", "windsp", "airtemp", 
                               "relhum", "pressure", "dailyrain", "par1", "par2",
                               "par3", "cdom", "mvolts", "co2-1", "co2-2", "temp1",
                               "cond1", "ph1", "ph1mV", "turb", "chl", "chlrfu", "bga1cell",
                               "bga1rfu", "ODOrel1", "ODOabs1", "temp2", "cond2", "ph2", 
                               "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3", 
                               "temp4", "temp5", "temp6", "temp7"))
testing <- read.csv("~/BPoundBuoy_2015_QC_Processed_Data.csv")
## Note that 'testing' follows the files sent by Julie on 15th Sep 2016, notes and data.
##    Note also that I used sheet C6 - this is the most up to date version at the minute

## Note email from Helen April 25, 2016: " 2015 both CO2 sensors failed after a few days due to 
##    a plistdip issue after spring deployment. Deep CO2 failed: 5/15/2015 10:50 (0-2000 ppm)
##    Shallow CO2 failed: 5/16/2015 3:10 (0-5000 ppm); 2 new sensors were put on
## They were put back out: 7/14/2015 10:50 But were probably not stabilized for a few hours
##    Deep sonde failed on:7/20/2015 13:10 (0-2000 ppm); The shallow sonde did not fail.
## Note: it was not cleaned all summer (0-2000ppm)"
## Note that ODO1 has no data for the functional CO2 perid =(
## FIXME: this sonde is still reading above 2000pm so is it still actually a 5000 sonde? ALSO!
##    The CO2 values are vastly above what the calculations imply... I find it hard to believe
##    them.... was something done wrong?

bdat <- transform(bdat, datetime = as.POSIXct(as.character(datetime), 
                                              format = "%m/%d/%Y %I:%M:%S %p"))
bdat <- transform(bdat, Hour = as.numeric(format(datetime, format = "%H")))
bdat <- transform(bdat, Month = as.numeric(format(datetime, format = "%m")))
bdat <- transform(bdat, Day = as.numeric(format(datetime, format = "%d")))
bdat <- transform(bdat, DOY = as.numeric(format(datetime, format = "%j")))
bdat <- transform(bdat, Time = as.numeric(format(datetime, "%H")) +
                            as.numeric(format(datetime, "%M"))/60)

testing <- transform(testing, DateTime = as.POSIXct(as.character(DateTime), 
                                                    format = "%Y-%m-%d %H:%M"))
testing <- transform(testing, Hour = as.numeric(format(DateTime, format = "%H")))
testing <- transform(testing, Month = as.numeric(format(DateTime, format = "%m")))
testing <- transform(testing, Day = as.numeric(format(DateTime, format = "%d")))
testing <- transform(testing, DOY = as.numeric(format(DateTime, format = "%j")))
testing <- transform(testing, Time = as.numeric(format(DateTime, "%H")) +
                       as.numeric(format(DateTime, "%M"))/60)

missing <- which(!testing$DateTime %in% bdat$datetime )
testing[missing,]# these are Missing in the master spreadsheet, must have been deleted in bdat original

missing2 <- which(!bdat$datetime %in% testing$DateTime)
bdat[missing2,] #these are the last dates of October that are not in the master, probs deleted due to 
#   probe being taken out?
testing <- testing[-missing,]
bdat <- bdat[-missing2,]

## ok now we have complete match of dates... let's check other vars!
plot(bdat$co2.1 ~ testing$CO2shallow, type = "l")
plot(bdat$ph1 ~ testing$pHShallow, type = "l")
plot(bdat$cdom ~ testing$CDOM, type = "l")

## yep that is good;now for flags that I might cause me to take more rows out or make NA
## 1. rest of October probe being taken out times
bdat[bdat$datetime >= as.POSIXct("06.10.2015 10:20", format="%d.%m.%Y %H:%M"),-1] <- NA
## 2. June times of initial deployment
bdat[bdat$datetime <= as.POSIXct("14.05.2015 16:30", format="%d.%m.%Y %H:%M"),-1] <- NA
## 3. Calibration periods
bdat[bdat$datetime >= as.POSIXct("11.06.2015 09:30", format="%d.%m.%Y %H:%M") &
       bdat$datetime <= as.POSIXct("11.06.2015 15:20", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat[bdat$datetime >= as.POSIXct("25.06.2015 10:40", format="%d.%m.%Y %H:%M") &
       bdat$datetime <= as.POSIXct("25.06.2015 14:10", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat[bdat$datetime >= as.POSIXct("09.07.2015 09:30", format="%d.%m.%Y %H:%M") &
       bdat$datetime <= as.POSIXct("09.07.2015 15:00", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat[bdat$datetime >= as.POSIXct("26.07.2015 13:10", format="%d.%m.%Y %H:%M") &
       bdat$datetime <= as.POSIXct("26.07.2015 14:20", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat[bdat$datetime >= as.POSIXct("05.08.2015 10:30", format="%d.%m.%Y %H:%M") &
       bdat$datetime <= as.POSIXct("05.08.2015 13:50", format="%d.%m.%Y %H:%M"),-1] <- NA
bdat[bdat$datetime >= as.POSIXct("03.09.2015 10:00", format="%d.%m.%Y %H:%M") &
       bdat$datetime <= as.POSIXct("03.09.2015 13:50", format="%d.%m.%Y %H:%M"),-1] <- NA


## let's look at when sun up and down
sundat <- sunrise.set(lat = 50.648016, long=-105.5072930, date = '2015/05/01', 
                      timezone = 'UTC+6', num.days=170)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

## merge sundat with BP data
bdat <- merge(bdat, sundat)

## isday to depend on sunup and down rather than rigid times
bdat$isDay <- ifelse(bdat$Time < bdat$DownTime & 
                               bdat$Time > bdat$UpTime, 
                             TRUE, FALSE)
bdat$TimeofDay <- ifelse(bdat$Time < bdat$DownTime & 
                                   bdat$Time > bdat$UpTime, 
                                 'Day', 'Night')

press <- read.csv("../data/bp-pressure2015.csv")
names(press)[which(names(press) == 'Pressure')] <- 'AirPressure'
press <- press[,-which(names(press) == 'Year')]
## get pressure in
bdat <- merge(bdat, press, by=c('Month', 'Day'))


## check the CO2 columns
with(bdat, plot(co2.1 ~ datetime))

## indeed june and july crazy; here's the density  
with(bdat, plot(density(co2.1), main='Density of 2015 uncorrected all-outliers-in pCO2'))
with(bdat, plot(density(co2.2)))

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

## make NA the dates when the shallow one was out of operation and put back in, 
##    and the few ones before June
outliers <- which(bdat$datetime < as.POSIXct("15/07/2015", format="%d/%m/%Y"))
bdat$co2.1[outliers] <- NA

## plot some diagnostics and look if pH/chl determining CO2
with(bdat, plot(co2.1 ~ Hour, col = ifelse(isDay, "red", "black")))

with(bdat, plot(co2.1 ~ ph1, col = ifelse(Month < 9, "red", "black")))

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
ml2015 <- read.csv("../data/maunaloa2.csv")
bdat$Year <- rep(2015)
bdat <- merge(bdat, ml2015[,c('Year', 'Month', 'pCO2interp')])
bdat <- transform(bdat, salcalc = salcalc(temp=temp1, cond=cond1, 
                                          dbar=(AirPressure + 0.8)/10))

bdat$dic <- rep(NA)
bdatflux <- with(bdat, gasExchangeUser(temp = temp1, cond = cond1, ph=ph1, 
                                       wind=windsp, alknotdic = FALSE, 
                                       salt = salcalc, kpa=AirPressure,
                                       pco2atm = pCO2interp, diccalc = TRUE))
bdatall <- cbind(bdat, bdatflux)

## save for later
saveRDS(bdatall, "../data/private/bpbuoy2015-mod.rds")

## plot data
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

ggplot(data=bdat, aes(y = co2.1, x = datetime, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  ylab("CO2 (ppm)") +
  xlab("Day of year 2015")

