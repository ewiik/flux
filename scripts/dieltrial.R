## looking at sonde data from Sept 2015
## diel patterns?
## units in source file: Y/M/D HH:MM:SS	C	uS	uS		%	mg/L	ppt
## can get some kind of sunrise sunset info here but not in tabular format for autodownload:
##    http://www.timeanddate.com/sun/canada/regina?month=9&year=2015

## load necessary packages
library("ggplot2")
library("scales")
library('StreamMetabolism')
library('mgcv')

## read in data
diel <- read.csv("../data/private/WS-9-9.csv")

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
diel <- transform(diel, Date.Time = as.POSIXct(as.character(Date.Time), 
                                               format = "%y/%m/%d %H:%M:%S"))
diel <- transform(diel, Day = as.numeric(format(Date.Time, format = "%d")))
diel <- transform(diel, Month = as.numeric(format(Date.Time, format = "%m")))
diel <- transform(diel, Time = as.numeric(format(Date.Time, "%H")) +
                            as.numeric(format(Date.Time, "%M"))/60)
diel <- transform(diel, DOY = as.numeric(format(Date.Time, "%j")))

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
diel <- transform(diel, isDay = ifelse(Hour < 21 & Hour > 9, TRUE, FALSE))

## let's look at when sun really up and down
sundat <- sunrise.set(lat = 50.429771, long=-104.608115, date = '2015/09/01', 
                      timezone = 'UTC+6', num.days=30)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

## merge
diel <- merge(diel, sundat)

## plot
with(diel, plot(Temp ~ Date.Time, col = ifelse(isDay, "red", "black")))

dielstack <- stack(diel[,2:8])
dielstack$Date.Time <- rep(diel$Date.Time)
dielstack$Day <- rep(diel$Day)

dielsplit <- with(dielstack, split(dielstack, list(ind)))

pH <- diel[,c('pH', 'Date.Time', 'isDay')]
names(pH)[1:2] <- c("y", "Date")
pH$panel <- rep("pH")
cond <- diel[,c('Cond', 'Date.Time', 'isDay')]
names(cond)[1:2] <- c("y", "Date")
cond$panel <- rep('COND')
temp <- diel[,c('Temp', 'Date.Time', 'isDay')]
names(temp)[1:2] <- c("y", "Date")
temp$panel <- rep("TEMP")
oxy <- diel[,c('ODOsat', 'Date.Time', 'isDay')]
names(oxy)[1:2] <- c("y", "Date")
oxy$panel <- rep("OXY")
oxymg <- diel[,c('ODO', 'Date.Time', 'isDay')]
names(oxymg)[1:2] <- c("y", "Date")
oxymg$panel <- rep("OXYmg")

super <- rbind(pH, cond, temp, oxy, oxymg)

p <- ggplot(data = super, mapping = aes(x = Date, y = y, color = isDay), size = 0.2) +
  scale_size_manual(values = c(0.2, 0.2)) +
  scale_color_manual(values = c("black", "darkgrey")) +
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = "grey", size = 0.5)) 
p <- p + facet_grid(panel~., scale="free")
p <- p + layer(data = pH, geom = c("point"), stat = "identity", position="identity")
p <- p + layer(data = cond, geom = c("point"), stat = "identity", position="identity")
p <- p + layer(data = oxy, geom = c("point"), stat = "identity", position="identity")
p <- p + layer(data = temp, geom = c("point"), stat = "identity", position="identity") 
p <- p + layer(data = oxymg, geom = c("point"), stat = "identity", position="identity") + 
  theme(legend.position = "top")
p <- p + geom_vline(data = sondeclean, aes(xintercept = as.numeric(Cleantimes), 
                                           shape = Cleaned),
                    show.legend = TRUE)
p

pdf("../data/private/dieltrial.pdf", width = 15)
p
dev.off()

## grab Sep 2015 mauna loa value
pco2atm <- ml[ml$Month == 9 & ml$Year == 2015, 'pCO2interp']

## need to create DIC from conductivity since we don't have instantaneous DIC or alk measures
## using the equation used by kerri for missing DIC values.
diel$dic <- 26.57 + 0.018 * diel$Cond  # (mg/L)
diel$dicumol <- diel$dic / 0.012 # --> uM

## run gasExchangeFlex, using mean wind (from kerri's means) and Wascana's altitude
CO2Fluxz <- with(diel, gasExchangeFlex(Temp, Cond, pH, wind = 2.8, kerri = FALSE, altnotkpa = TRUE,
                                       salt = Sal, dic = dicumol, alt = 570.5, pco2atm = pco2atm))
diel <- transform(diel, CO2Flux = CO2Fluxz$fluxenh)
diel <- transform(diel, pCO2 = CO2Fluxz$pco2)

## update isDay to the sun up thing
diel <- transform(diel, isDay = ifelse(Time < DownTime & Time > UpTime, TRUE, FALSE))


## save object for other purposes
saveRDS(diel, '../data/private/WS-9-9-mod.rds')

with(diel, plot(CO2Flux ~ Date.Time, col = ifelse(isDay, 'red', 'black')))

ggplot(data=diel, aes(y = pCO2, x = Time, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  scale_color_identity() + # means it understands red and black in ifelse
  scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  #geom_vline(xintercept = as.numeric(datesDOY)) +
  ylab("CO2 (ppm)") +
  xlab("Hour")

ggplot(data=diel, aes(y = Temp, x = Date.Time, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  #geom_vline(xintercept = as.numeric(datesDOY)) +
  ylab("Temperature") +
  xlab("Hour")

ggplot(bdat2014fullf, aes(isDay, co2corr, group=Hour)) +
  geom_boxplot()

## gamming the cyclic stuff in the data
m1 <- gam(pH ~ ti(Time, bs = "cc") + ti(DOY), data = diel)

m2 <- gam(pH ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")), 
                                                  data = diel)
anova(m1, m2, test = "LRT")

plot(m2, scheme=2, ylab='mean-0 pH', main='mean-0 pH') # latter works for 3D plot
#   former for the 2D plot
plot(acf(resid(m2)))

pdf("../data/private/diel-wwgam.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(m2, pages = 4, scheme = 2)
par(op)
dev.off()

## predict for periods of varying interactions.
N <- 200
simDOY <- c(254:256, 262:264, 270:272)
DOYgroup <- factor(rep(1:3, each=3))
DOYgroup <- factor(rep(c('254:256','262:264','270:272'), each=3))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(diel$Time, na.rm=TRUE),max(diel$Time, na.rm=TRUE),
                                  length = N), times=reptimes),
           `DOY` = rep(simDOY, each = N))
m2pred <- predict(m2, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, m2pred, DOYgroups)
names(predicted)[which(names(predicted)=='m2pred')] <- 'pH'
ggplot(predicted, aes(x = Time, y = pH, group= DOY, col=DOYgroups)) +
  theme_bw() +
  geom_line() +
  scale_colour_discrete(name="Day of Year") +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab('pH')

