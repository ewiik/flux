## load necessary packages
library("ggplot2")
library('StreamMetabolism')
library('mgcv')
library('viridis')
library('reshape2')
## Separate script for raw data in 2016, received from Helen by email 29.3.2017
## load necessary functions
source('../functions/gasExchangeUser.R')

if (!file.exists('../data/private/salcalc.R')) {
  stop("Get salinity calculation function from Emma")}
source('../data/private/salcalc.R')

## load in supplementary data
if (!file.exists('../data/private/regvars.rds')) {
  source('../scripts/regression_routines.R')
}
regvars <- readRDS('../data/private/regvars.rds')
if (!file.exists('../data/private/params-flux.rds')) {
  source("../scripts/co2_scenarios.R")
}
params <- readRDS('../data/private/params-flux.rds')

## mauna loa stuff online has been updated from my getmaunaloa script
if (!file.exists("../data/maunaloa2.csv")) {
  download.file("ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt", 
                destfile = "../data/maunaloa2.csv") 
  ml <- read.csv("../data/maunaloa2.csv", skip = 72, sep = "", encoding = "latin1", header = FALSE,
                 col.names = c("Year", "Month", "DecimalDate", "pCO2ave", "pCO2interp", 
                               "trend", "hashdays")) 
  write.csv(ml, "../data/maunaloa2.csv", row.names = FALSE)
}
ml <- read.csv("../data/maunaloa2.csv")

## raw data
bdat2016 <- read.csv("../data/private/BP2016RawData - Copy.csv", skip = 4,
                         col.names = c('datetime',"batvolt", "winddir", "windsp", "airtemp",
                                       "relhum", "pressure", "dailyrain", "par1", "par2",
                                       "par3", "cdom", "mvolts", "co2-1", "co2-2", "temp1",
                                       "spcond1", "ph1", "ph1mV", "turb", "chl", "chlrfu", "bga1cell",
                                       "bga1rfu", "ODOrel1", "ODOabs1", "temp2", "spcond2", "ph2",
                                       "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3",
                                       "temp4", "temp5", "temp6", "temp7"))
## make time understandable
bdat2016 <- transform(bdat2016, datetime = as.POSIXct(as.character(datetime),
                                                      format = "%m/%d/%Y %r",  tz="Canada/Saskatchewan")) 
bdat2016 <- transform(bdat2016, Hour = as.numeric(format(datetime, format = "%H")))
bdat2016 <- transform(bdat2016, Month = as.numeric(format(datetime, format = "%m")))
bdat2016 <- transform(bdat2016, Day = as.numeric(format(datetime, format = "%d")))
bdat2016 <- transform(bdat2016, DOY = as.numeric(format(datetime, format = "%j")))
bdat2016 <- transform(bdat2016, Time = as.numeric(format(datetime, "%H")) +
                        as.numeric(format(datetime, "%M"))/60)

## make factor numerics intno numerics
facs <- c("airtemp","winddir","windsp","relhum", "pressure", "dailyrain")
bdat2016[,facs] <- apply(bdat2016[,facs], 2, as.character)

makenas <- apply(bdat2016[,facs], 2, function(x) grep("-Invalid", x))
faclist <- as.list(bdat2016[,facs])
for (i in 1:length(faclist)) {
  faclist[[i]][makenas[[i]]] <- NA
  bdat2016[,facs[i]] <- faclist[[i]]
}
bdat2016[,facs] <- apply(bdat2016[,facs], 2, as.numeric)


## make -100000 and the likes into NA
meters <- c("batvolt" , "winddir", "windsp", "airtemp", "relhum", "pressure",  "dailyrain", "par1",
            "par2", "par3", "cdom", "mvolts", "co2.1", "co2.2", "temp1", "spcond1",  "ph1", "ph1mV",
            "turb", "chl", "chlrfu", "bga1cell", "bga1rfu",  "ODOrel1", "ODOabs1", "temp2", "spcond2",
            "ph2", "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3", "temp4", "temp5", 
            "temp6", "temp7" )
makenas <- apply(bdat2016[,meters], 2, function(x) which(x<= -5000))
meterlist <- as.list(bdat2016[,meters])
for (i in 1:length(meterlist)) {
  meterlist[[i]][makenas[[i]]] <- NA
  bdat2016[,meters[i]] <- meterlist[[i]]
}
mapply(function(x,y) plot(density(x, na.rm=TRUE), main=y), bdat2016[,meters],meters)
# outliers in temp7,temp6,temp5,temp4,temp3 (>30), ODOabs, bga2rfu (<0), ph2, spcon2 (<700), temp2, 
# ODOabs1, bga1rfu, potentially turb, ph1, spcond1, temp1, co2.2 (peak 0 at bimodal.. <0?), and aroun 2000...
#   reliable?, par3 heavy right tail, par2 -ve values, pressure (>100), 

## most weirdness occurs early May.... remove these dates altogether
bdat2016 <- bdat2016[-which(bdat2016$DOY < 140),]
# outliers remaining: ph2, 

## remaining outliers:
bdat2016$ph2[bdat2016$ph2 < 7.9] <- NA
bdat2016$spcond1[bdat2016$spcond1 < 740] <- NA
bdat2016$spcond2[bdat2016$spcond2 < 600] <- NA #basically a few near-0s
bdat2016$pressure[bdat2016$pressure > 200] <- NA 
bdat2016$co2.1[bdat2016$co2.1 >= 1900] <- NA 


## add pressure in kpa (original in inHg)
bdat2016$pressureKPA <- bdat2016$pressure * 3.3864
bdat2016$pressureKPAwater <- bdat2016$pressureKPA + .08

## calculate salinity
bdat2016 <- transform(bdat2016, salcalc = salcalc(temp=temp1, cond=spcond1,
                                                          dbar=pressureKPAwater/10))


## calculate day lengths etc
## let's look at when sun up and down
sundat <- sunrise.set(lat = 50.648016, long=-105.5072930, date = '2014/06/01',
                      timezone = "Canada/Saskatchewan", num.days=90)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

## merge sundat with BP data
bdat2016 <- merge(bdat2016, sundat)

## make isday index based on sunup and down
bdat2016$isDay <- ifelse(bdat2016$Time < bdat2016$DownTime &
                               bdat2016$Time > bdat2016$UpTime,
                             TRUE, FALSE)
bdat2016$TimeofDay <- ifelse(bdat2016$Time < bdat2016$DownTime &
                                   bdat2016$Time > bdat2016$UpTime,
                                 'Day', 'Night')
## merge with maunaloa
bdat2016$Year <- rep(2016)
bdat2016 <- merge(bdat2016, ml[,c("Year", "Month", "pCO2interp")])

## create correct units for spcond and windsp (mph), and correct co2
bdat2016$windms <- bdat2016$windsp * 0.44704
bdat2016$cond1us <- bdat2016$spcond1 * (1 + 0.02 * (bdat2016$temp1 - 25)) 
bdat2016$cond2us <- bdat2016$spcond2 * (1 + 0.02 * (bdat2016$temp1 - 25)) 

corrpress <- bdat2016$co2.1*((1013 - bdat2016$pressureKPA*10) * 0.0015)
corrtemp <- bdat2016$co2.1*((25 - bdat2016$temp1) * 0.003)
corrairpress <- -bdat2016$co2.1*((1014.14 - 1013.25) * 0.0015)
# 89cm pressure calculated at http://www.calctool.org/CALC/other/games/depth_press
#   is where the CO2 sensor is.

bdat2016$co2.1corr <- bdat2016$co2.1 + corrpress + corrtemp + corrairpress

## calculate flux
bdat2016flux <- with(bdat2016, gasExchangeUser(temp = temp1, cond = cond1us, ph=ph1,
                                                       wind=windms, alknotdic = FALSE,
                                                       salt = salcalc, kpa=pressureKPA,
                                                       pco2atm = pCO2interp, diccalc = TRUE))
bdat2016all <- cbind(bdat2016, bdat2016flux)
names(bdat2016all)[which(names(bdat2016all)=="pco2")] <- "pco2calc"


## plot some diagnostics
with(bdat2016all, plot(co2.1corr ~ datetime))
with(bdat2016all[order(bdat2016all$datetime),], lines(pco2calc ~ datetime, col='green'))
