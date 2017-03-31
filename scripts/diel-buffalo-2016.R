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
testing <- mapply(function(x,y) x[y] <- NA, faclist, makenas)
## FIXME: this works faclist[[1]][makenas[[1]]] but not the mapply!
## FIXME: proceed to make stuff below work
bdat2016[makenas,facs] <- NA

bdat2016[,facs] <- apply(bdat2016[,facs], 2, as.numeric)


## calculate salinity
bdat2016 <- transform(bdat2016, salcalc = salcalc(temp=temp1, cond=cond1,
                                                          dbar=pressureKPAwater/10))
bdat2016flux <- with(bdat2016, gasExchangeUser(temp = temp1, cond = cond1, ph=ph1,
                                                       wind=windsp, alknotdic = FALSE,
                                                       salt = salcalc, kpa=pressureKPA,
                                                       pco2atm = maunaloa, diccalc = TRUE))
bdat2016all <- cbind(bdat2016, bdat2016flux)
bdat2016corr <- subset(bdat2016, select=c('datetime', 'co2corr'))
bdat2016full <- merge(bdat2016all, bdat2016corr, all.x = TRUE)


## take away nonsense values


