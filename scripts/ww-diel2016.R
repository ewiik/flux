## dealing with Wascana 2016 data

## load necessary packages
library("ggplot2")
library("scales")
library('StreamMetabolism')
library('mgcv')
library('viridis')
library('extrafont')

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
diel <- read.csv("../data/private/Wascana_2016_complete.csv", skip = 2, header=FALSE)
names(diel) <- c("Date.Time","temp","condus", "spcondus","ph","odorel","odomgl","salppt")

# get Kyle's station data (see email July 5, 2015)
if (!file.exists("../data/hodder20164.csv")) {
  urls <- c(paste0("http://uregina.ca/~hodder2k/data/2016_0",4:9,".txt"), 
            "http://uregina.ca/~hodder2k/data/2016_10.txt")
  mapply(FUN=download.file, urls, destfile = paste0("../data/hodder2016",4:10,".csv"))
}
if (!file.exists("../data/hodder2016.rds")) {
  readfiles = paste0("../data/hodder2016",4:10,".csv")
  kyles <- lapply(readfiles, read.csv, header = FALSE)
  kyles <- do.call(rbind, kyles)
  names(kyles) <- c("date", "time","airtemp","relhum","dewpoint","windspeedkmh","highwindkmh","avwindbear","rainmmh",
                   "cumrainmm","hpa","totrainmm","rictemp","ricrelhum","windgustkmh","windchill","heatind",
                   "NA1","NA2","NA3","NA4","apptemp","maxrad","cumsunhours","windbear","NA5","totrainmidnightmm")
  kyles$datetime <- paste(kyles$date, kyles$time)
  kyles <- transform(kyles, date=as.POSIXlt(as.character(date), format="%d/%m/%y", tz="Canada/Saskatchewan"))
  kyles <- transform(kyles, datetime=as.POSIXlt(as.character(datetime), format="%d/%m/%y %H:%M", 
                                              tz="Canada/Saskatchewan"))
  saveRDS(kyles, "../data/hodder2016.rds")
}
kyles <- readRDS("../data/hodder2016.rds")

kyles$windms <- kyles$windspeedkmh * (1000/60/60)
kyles$kpa <- kyles$hpa/10

## load functions
source("../functions/gasExchangeFlex.R")

## mauna loa stuff online has been updated from my getmaunaloa script, and since I need
##    2016....
if (!file.exists("../data/maunaloa2.csv")) {
  download.file("ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt", 
                destfile = "../data/maunaloa2.csv") 
  ml <- read.csv("../data/maunaloa2.csv", skip = 72, sep = "", encoding = "latin1", header = FALSE,
                 col.names = c("Year", "Month", "DecimalDate", "pCO2ave", "pCO2interp", 
                               "trend", "hashdays")) 
  write.csv(ml, "../data/maunaloa2.csv", row.names = FALSE)
}
ml <- read.csv("../data/maunaloa2.csv")

## transform date, remove data points prior to instrument stabilisation
diel <- transform(diel, Date.Time = as.POSIXlt(as.character(Date.Time), 
                                               format = "%y/%m/%d %H:%M:%S",
                                               tz="Canada/Saskatchewan"))
diel$Date.Time <- round(diel$Date.Time, units = "mins")
diel <- transform(diel, Day = as.numeric(format(Date.Time, format = "%d")))
diel <- transform(diel, Month = as.numeric(format(Date.Time, format = "%m")))
diel <- transform(diel, Time = as.numeric(format(Date.Time, "%H")) +
                    as.numeric(format(Date.Time, "%M"))/60)
diel <- transform(diel, DOY = as.numeric(format(Date.Time, "%j")))
diel <- transform(diel, Hour = as.numeric(format(Date.Time, format = "%H")))

## let's look at when sun really up and down
sundat <- sunrise.set(lat = 50.429771, long=-104.608115, date = '2016/04/01', 
                      timezone = 'UTC+6', num.days=250)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

## merge diel with sundat and kyle (all diel timestamps are contained in kyle)
## ps merge didn't seem to work with posixlt cause they lists, so made into char (tested
##    which(!test3$Date.Time == test3$datetime) to check all was ok upon merging by char
## Note that lt rather than ct time object initially chosen cause ct made seconds unmatchable for diel
diel <- merge(diel, sundat)

take <- which(kyles$datetime %in% diel$Date.Time)
kyles$datchar <- as.character(kyles$datetime)
diel$datchar <- as.character(diel$Date.Time)

diel <- merge(diel, kyles, by="datchar")
takeout <- which(names(diel) %in% c("datchar", "Date.Time"))
diel <- diel[,-takeout]
diel$datetime <- as.POSIXct(diel$datetime)

## create day or night index
diel <- transform(diel, isDay = ifelse(Time < DownTime & Time > UpTime, TRUE, FALSE))

## reorder to time
diel <- diel[order(diel$datetime),]

## merge with maunaloa
diel$Year <- rep(2016)
diel <- merge(diel, ml[,c('Year','Month','pCO2interp')])
names(diel)[grep("interp", names(diel))] <- "pCO2ml"

## need to create DIC from conductivity since we don't have instantaneous DIC or alk measures
## using the equation used by kerri for missing DIC values.
diel$dic <- 26.57 + 0.018 * diel$condus  # (mg/L)
diel$dicumol <- diel$dic / 0.012 # --> uM

## take out nonsense values
diel <- diel[-which(diel$ph < 8.2),] # removes still-in-ph4buffer times, and sonde blips
diel$condus[which(diel$condus < 800)] <- NA # some malfunction but other probes ok
diel$salppt[which(diel$condus < 800)] <- NA
diel$salppt[which(diel$salppt < 0.5)] <- NA

## FIXME: take out sonde removal readings once date from D
#diel <- diel[-which()]

## run gasExchangeFlex, using Kyle's recorded wind and pressure
CO2Fluxz <- with(diel, gasExchangeFlex(temp, condus, ph, wind = windms, kerri = FALSE, altnotkpa = FALSE,
                                       salt = salppt, dic = dicumol, kpa = kpa, pco2atm = pCO2ml))
diel <- transform(diel, CO2Flux = CO2Fluxz$fluxenh)
diel <- transform(diel, pCO2 = CO2Fluxz$pco2)

## save rds for further examination and modeling
saveRDS(diel, '../data/private/diel-W-2016.rds')
