## looking at sonde data from Sept 2015
## units in source file: Y/M/D HH:MM:SS	C	uS	uS		%	mg/L	ppt

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
diel <- read.csv("../data/private/WS-9-9.csv")
# get Kyle's station data from September 2015 (see email July 5, 2015)
if (!file.exists("../data/hodderSep2015.csv")) {
  download.file("http://uregina.ca/~hodder2k/data/2015_09.txt", destfile = "../data/hodderSep2015.csv")
}
if (!file.exists("../data/hodderSep2015.rds")) {
kyle <- read.csv("../data/hodderSep2015.csv", header = FALSE)
names(kyle) <- c("date", "time","airtemp","relhum","dewpoint","windspeedkmh","highwindkmh","avwindbear","rainmmh",
                 "cumrainmm","hpa","totrainmm","rictemp","ricrelhum","windgustkmh","windchill","heatind",
                 "NA1","NA2","NA3","NA4","apptemp","maxrad","cumsunhours","windbear","NA5","totrainmidnightmm")
kyle$datetime <- paste(kyle$date, kyle$time)
kyle <- transform(kyle, date=as.POSIXlt(as.character(date), format="%d/%m/%y", tz="Canada/Saskatchewan"))
kyle <- transform(kyle, datetime=as.POSIXlt(as.character(datetime), format="%d/%m/%y %H:%M", 
                                            tz="Canada/Saskatchewan"))
saveRDS(kyle, "../data/hodderSep2015.rds")
# http://uregina.ca/~hodder2k/archive.htm
}
kyle <- readRDS("../data/hodderSep2015.rds")
kyle$windms <- kyle$windspeedkmh * (1000/60/60)
kyle$kpa <- kyle$hpa/10

## load functions
source("../functions/gasExchangeFlex.R")

## mauna loa stuff online has been updated from my getmaunaloa script, and since I need
##    2015....
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

takeout <- which(diel$Cond < 100) 
tocheck <- diel[takeout,] # we have some abnormalities 14th of Sept!! between 13:25 & 13:55
#   when the sonde was taken out and cleaned etc. 
diel <- diel[-takeout,]
tocheck <- transform(tocheck, Day = as.numeric(format(Date.Time, format = "%d")))
tocheck <- transform(tocheck, Month = as.numeric(format(Date.Time, format = "%m")))
cleaningtime <- which(tocheck$Day == 14 & tocheck$Month == 9)
sondeclean <- data.frame(Cleantimes = c(tocheck[cleaningtime[1],1], 
                                        tocheck[cleaningtime[length(cleaningtime)],1]), Cleaned = "")

## let's look at when sun really up and down
sundat <- sunrise.set(lat = 50.429771, long=-104.608115, date = '2015/09/01', 
                      timezone = 'UTC+6', num.days=30)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

## merge diel with sundat and kyle (all diel timestamps are contained in kyle)
## ps merge didn't seem to work with posixlt cause they lists, so made into char (tested
##    which(!test3$Date.Time == test3$datetime) to check all was ok upon merging by char
## Note that lt rather than ct time object initially chosen cause ct made seconds unmatchable for diel
diel <- merge(diel, sundat)

take <- which(kyle$datetime %in% diel$Date.Time)
kyle$datchar <- as.character(kyle$datetime)
diel$datchar <- as.character(diel$Date.Time)

diel <- merge(diel, kyle, by="datchar")
takeout <- which(names(diel) %in% c("datchar", "Date.Time"))
diel <- diel[,-takeout]
diel$datetime <- as.POSIXct(diel$datetime)

## create day or night index
diel <- transform(diel, isDay = ifelse(Time < DownTime & Time > UpTime, TRUE, FALSE))

## reorder to time
diel <- diel[order(diel$datetime),]

## grab Sep 2015 mauna loa value
pco2atm <- ml[ml$Month == 9 & ml$Year == 2015, 'pCO2interp']

## need to create DIC from conductivity since we don't have instantaneous DIC or alk measures
## using the equation used by kerri for missing DIC values.
diel$dic <- 26.57 + 0.018 * diel$Cond  # (mg/L)
diel$dicumol <- diel$dic / 0.012 # --> uM

## run gasExchangeFlex, using Kyle's recorded wind and pressure
CO2Fluxz <- with(diel, gasExchangeFlex(Temp, Cond, pH, wind = windms, kerri = FALSE, altnotkpa = FALSE,
                                       salt = Sal, dic = dicumol, kpa = kpa, pco2atm = pco2atm))
diel <- transform(diel, CO2Flux = CO2Fluxz$fluxenh)
diel <- transform(diel, pCO2 = CO2Fluxz$pco2)

## take away times where sonde was pulled out of the water
diel <- diel[-which(diel$datetime >= 
                      as.POSIXct("2015-09-29 14:05", format = "%Y-%m-%d %H:%M")),]

## save object for other purposes
saveRDS(diel, '../data/private/WS-9-9-mod.rds')

## plot
dielstack <- stack(diel[,2:8])
dielstack$datetime <- rep(diel$datetime)
dielstack$Day <- rep(diel$Day)

dielsplit <- with(dielstack, split(dielstack, list(ind)))

pH <- diel[,c('pH', 'datetime', 'isDay')]
names(pH)[1:2] <- c("y", "Date")
pH$panel <- rep("pH")
cond <- diel[,c('Cond', 'datetime', 'isDay')]
names(cond)[1:2] <- c("y", "Date")
cond$panel <- rep('COND')
temp <- diel[,c('Temp', 'datetime', 'isDay')]
names(temp)[1:2] <- c("y", "Date")
temp$panel <- rep("TEMP")
oxy <- diel[,c('ODOsat', 'datetime', 'isDay')]
names(oxy)[1:2] <- c("y", "Date")
oxy$panel <- rep("OXY")
oxymg <- diel[,c('ODO', 'datetime', 'isDay')]
names(oxymg)[1:2] <- c("y", "Date")
oxymg$panel <- rep("OXYmg")

super <- rbind(pH, cond, temp, oxy, oxymg)

p <- ggplot(data = super, mapping = aes(x = Date, y = y, color = isDay), size=0.1) +
  scale_color_viridis(discrete = TRUE, alpha=0.4, begin=0.3, end=0.8, name="", breaks=
                        c(TRUE, FALSE), labels=c("Day", "Night")) +
  papertheme + 
  theme(axis.title.y=element_blank())  
p <- p + facet_grid(panel~., scale="free")
p <- p + layer(data = pH, geom = c("point"), stat = "identity", position="identity")
p <- p + layer(data = cond, geom = c("point"), stat = "identity", position="identity")
p <- p + layer(data = oxy, geom = c("point"), stat = "identity", position="identity")
p <- p + layer(data = temp, geom = c("point"), stat = "identity", position="identity") 
p <- p + layer(data = oxymg, geom = c("point"), stat = "identity", position="identity") 
p <- p + geom_vline(data = sondeclean, aes(xintercept = as.numeric(Cleantimes), 
                                           shape = Cleaned),show.legend = FALSE)
p

pdf("../docs/private/dieltrial.pdf", width = 15)
p
dev.off()

ggplot(data=diel, aes(y = pCO2, x = datetime, col = isDay)) + 
  scale_color_viridis(discrete = TRUE, alpha=0.4, begin=0.3, end=0.8, name="", breaks=
                        c(TRUE, FALSE), labels=c("Day", "Night")) +
  papertheme + 
  geom_point() +
  ylab("CO2 (ppm)") +
  xlab("Hour")

ggplot(data =diel, aes(y=Temp, x=datetime, col = isDay)) +
  papertheme +
  geom_point(size=0.4) +
  scale_color_viridis(alpha=0.4, begin=0.3, end=0.8, discrete = TRUE, name="", breaks =
                        c(TRUE, FALSE), labels=c("Day", "Night")) +
  ylab("Temperature") + xlab("Date")

