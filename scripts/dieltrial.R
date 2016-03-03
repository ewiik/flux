## looking at sonde data from Sept 2015
## diel patterns?
## units in source file: Y/M/D HH:MM:SS	C	uS	uS		%	mg/L	ppt
## can get some kind of sunrise sunset info here but not in tabular format for autodownload:
##    http://www.timeanddate.com/sun/canada/regina?month=9&year=2015

## load necessary packages
library("ggplot2")
library("scales")

## read in data, transform date, remove data points prior to instrument stabilisation
diel <- read.csv("../data/private/WS-9-9.csv")
diel <- transform(diel, Date.Time = as.POSIXct(as.character(Date.Time), 
                                               format = "%y/%m/%d %H:%M:%S"))
takeout <- which(diel$Cond < 100) 
tocheck <- diel[takeout,] # we have some abnormalities 14th of Sept!! between 13:25 & 13:55
#   when the sonde was taken out and cleaned etc. 
diel <- diel[-takeout,]
tocheck <- transform(tocheck, Day = as.numeric(format(Date.Time, format = "%d")))
tocheck <- transform(tocheck, Month = as.numeric(format(Date.Time, format = "%m")))
cleaningtime <- which(tocheck$Day == 14 & tocheck$Month == 9)
sondeclean <- data.frame(Cleantimes = c(tocheck[cleaningtime[1],1], 
                         tocheck[cleaningtime[length(cleaningtime)],1]), Cleaned = "")
## create day or night index
diel <- transform(diel, Hour = as.numeric(format(Date.Time, format = "%H")))
diel <- transform(diel, Day = ifelse(Hour < 21 & Hour > 9, TRUE, FALSE))

## plot
with(diel, plot(Temp ~ Date.Time, col = ifelse(Day, "red", "black")))

dielstack <- stack(diel[,2:8])
dielstack$Date.Time <- rep(diel$Date.Time)
dielstack$Day <- rep(diel$Day)

dielsplit <- with(dielstack, split(dielstack, list(ind)))

pH <- diel[,c('pH', 'Date.Time', 'Day')]
names(pH)[1:2] <- c("y", "Date")
pH$panel <- rep("pH")
cond <- diel[,c('Cond', 'Date.Time', 'Day')]
names(cond)[1:2] <- c("y", "Date")
cond$panel <- rep('COND')
temp <- diel[,c('Temp', 'Date.Time', 'Day')]
names(temp)[1:2] <- c("y", "Date")
temp$panel <- rep("TEMP")
oxy <- diel[,c('ODOsat', 'Date.Time', 'Day')]
names(oxy)[1:2] <- c("y", "Date")
oxy$panel <- rep("OXY")
oxymg <- diel[,c('ODO', 'Date.Time', 'Day')]
names(oxymg)[1:2] <- c("y", "Date")
oxymg$panel <- rep("OXYmg")

super <- rbind(pH, cond, temp, oxy, oxymg)

p <- ggplot(data = super, mapping = aes(x = Date, y = y, color = Day), size = 0.2) +
              scale_size_manual(values = c(0.2, 0.2)) +
              scale_color_manual(values = c("black", "darkgrey")) +
              theme_bw() + 
              theme(panel.grid.major = element_line(colour = "grey", size = 0.5)) 
p <- p + facet_grid(panel~., scale="free")
p <- p + layer(data = pH, geom = c("point"), stat = "identity")
p <- p + layer(data = cond, geom = c("point"), stat = "identity")
p <- p + layer(data = oxy, geom = c("point"), stat = "identity")
p <- p + layer(data = temp, geom = c("point"), stat = "identity") 
p <- p + layer(data = oxymg, geom = c("point"), stat = "identity") + theme(legend.position = "top")
p <- p + geom_vline(data = sondeclean, aes(xintercept = as.numeric(Cleantimes), 
                                           shape = Cleaned),
                                           show_guide = TRUE)
p

pdf("../data/private/dieltrial.pdf", width = 15)
p
dev.off()

## let's chech how CO2 flux gets modelled by this.
source("../functions/gasExchangeFlex.R")

if (!file.exists("../data/maunaloa.csv")) {
  source("../functions/getmaunaloa.R")
}

## mauna loa stuff online has been updated from my getmaunaloa script, and since I need
##    2015....
download.file("ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt", 
              destfile = "../data/maunaloa2.csv")
ml <- read.csv("../data/maunaloa2.csv", skip = 72, sep = "", encoding = "latin1", header = FALSE,
               col.names = c("Year", "Month", "DecimalDate", "pCO2ave", "pCO2interp", 
                             "trend", "hashdays"))
write.csv(ml, "../data/maunaloa2.csv", row.names = FALSE)

## grab Sep 2015 value
pco2atm <- ml[ml$Month == 9 & ml$Year == 2015, 'pCO2interp']

## need to create DIC from conductivity since we don't have instantaneous DIC or alk measures
## using the equation used by kerri for missing DIC values.
dic <- 26.57 + 0.018 * diel$Cond  # (mg/L)
dic <- dic / 0.012 # --> uM

## run gasExchangeFlex, using mean wind (from kerri's means) and Wascana's altitude
source("../functions/gasExchangeFlex.R")
CO2Flux <- with(diel, gasExchangeFlex(Temp, Cond, pH, wind = 2.8, kerri = FALSE, altnotkpa = TRUE,
                                                  salt = Sal, dic = dic, alt = 570.5, pco2atm = pco2atm))
diel <- transform(diel, CO2Flux = CO2Flux$fluxenh)

with(diel, plot(CO2Flux ~ Date.Time, col = ifelse(Day, "red", "black")))
