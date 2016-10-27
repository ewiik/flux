## let's look at background CO2 levels in Saskatchewan/prairies: data received from 
##    Environment Canada **personally** -- see email 13.10.2016 (Esther and Bratt
##    did not exist at the links to public repos)
## NB!! For 2014, where we would like data to support diel Buffalo Pound story, there are
##    no data for Trout June-July; Esther Aug-Sep; Bratt Jul-Sep --> borrow data from other 
##    years?
## Bratt just 20km south of Regina (51°12'N, 104°42'W); Esther in the prairies just on the
##    Alberta side (51°40'N, 110°12'W); Trout is in the boreal zone (54°21'N, 104°59'W); 
##    see https://www.ec.gc.ca/mges-ghgm/default.asp?lang=En&n=87F5053F-1

## load packages
library("ggplot2")
library("StreamMetabolism")
library("mgcv")
library("viridis")
library("reshape2")
library("extrafont")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=18, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
if (!file.exists("../data/private/Esther-CO2-Hourly.DAT")) {
  stop("get CO2 data files from Emma or EC")}
esther <- read.delim("../data/private/Esther-CO2-Hourly.DAT", sep = ",")
bratt <- read.delim("../data/private/Bratts_Lake-CO2-Hourly.DAT", sep = ",")
trout <- read.delim("../data/private/East_Trout_Lake-CO2-Hourly.DAT", sep = ",")

if (!file.exists("../data/dmet.rds")) {
source("../scripts/gasexchange_pressuredata.R")}
temp <- readRDS("../data/dmet.rds")
reg <- subset(temp, StationID==3002 | StationID==51441) # subsetted with Bratt in mind
reg$Time <- as.integer(reg$Time)

## put them all together
bratt$Site <- "Bratt"
esther$Site <- "Esther"
trout$Site <- "Trout"

bet <- rbind(bratt, esther, trout)

## make a Date and Month column
bet$datetime <- with(bet, paste(Year, DayofYear, HourEnding, sep="-"))
bet$datetime <- as.POSIXct(bet$datetime, format = "%Y-%j-%H", tz="Canada/Saskatchewan")
bet$datetime <- bet$datetime - 30*60
bet$Month <- as.numeric(format(bet$datetime, "%m"))
bet$Day <- as.numeric(format(bet$datetime, "%d"))
bet$Date <- format(bet$datetime, format="%m-%d")

## create Day or not column
## let's look at when sun up and down
sundat <- sunrise.set(lat = 50.159721, long= -104.646614, date = '2010/01/01', 
                      timezone = "Canada/Saskatchewan", num.days=2555)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, Year = as.numeric(format(sunrise, format = "%Y")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

bet <- merge(bet, reg, by.x=c("Year", "Month", "Day", "HourEnding"), 
             by.y=c("Year", "Month", "Day","Time"))
names(bet)[grep("degC", names(bet))] <- c("temp", "dewtemp")

bet <- merge(bet, sundat)
bet$TOD <- ifelse(bet$HourEnding-0.5 < bet$DownTime & 
                    bet$HourEnding-0.5 > bet$UpTime, 
                  "Day", "Night")

## plot them and look
brattplot <- ggplot(data=bet[bet$Year>=2009 & bet$Site == "Bratt",], aes(y=Mean, x=DayofYear, col=TOD)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_wrap("Year", ncol=3) +
  scale_color_viridis(discrete=TRUE, option = "magma", direction = -1,
                      begin=0, end=0.7)  +
  #scale_color_viridis(name="Hour") + for when plotting by Hour
  theme(legend.key.width=unit(2,"cm")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)), legend.title=element_blank()) +
  ylab(expression(paste(italic('p')*'CO'[2]))) + xlab("Day of Year (200 ca mid-July)")

ggplot(data=bet[bet$Year>=2009 & bet$Site == "Esther",], aes(y=Mean, x=DayofYear, col=TOD)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_wrap("Year", ncol=4) +
  scale_color_viridis(discrete=TRUE, option = "magma", direction = -1,
                      begin=0, end=0.7)  +
  #scale_color_viridis(name="Hour") + for when plotting by Hour
  theme(legend.key.width=unit(2,"cm"), legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylab("CO2 (ppm)") + xlab("Day of Year (200 mid-July)")

ggplot(data=bet, aes(y=Mean, x=DayofYear, col=TOD)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_grid(Year ~ Site) +
  scale_color_viridis(discrete=TRUE, option = "magma", direction = -1,
                      begin=0, end=0.7)  +
  #scale_color_viridis(name="Hour") + for when plotting by Hour
  theme(legend.key.width=unit(2,"cm"), legend.title=element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylab("CO2 (ppm)") + xlab("Day of Year (200 mid-July)")

ggplot(data=bet[bet$Year>=2010& bet$Site =="Esther",], aes(y=Mean, x=DayofYear, col=HourEnding)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_wrap("Year", ncol=4) +
  scale_color_viridis() +
  theme(legend.key.width=unit(2,"cm")) 

ggplot(data=bet[bet$Year>=2010& bet$Site == "Trout",], aes(y=Mean, x=DayofYear, col=HourEnding)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_wrap("Year") +
  scale_color_viridis(name="Hour")  +
  theme(legend.key.width=unit(2,"cm"))

## model daily fluctuations in summer borrowing strength from other years?
brattmodnull <- gam(Mean ~ ti(DayofYear, bs="cc", k=80) + s(Year, k=5) + ti(HourEnding, bs="cc",k=10) +
                  ti(DayofYear, HourEnding, bs = c("cc","cc")), 
                data=bet[bet$Site == "Bratt",], family=scat(),
                method="REML", select=TRUE)
brattmod <- gam(Mean ~ s(Year, k=5) + te(DayofYear, HourEnding, bs = c("cc","cc"), k=c(40,15)), 
                data=bet[bet$Site == "Bratt",], family=scat(),
                method="REML", select=TRUE)

tempmod <- gam(temp ~ s(Year, k=5) + te(DayofYear, HourEnding, bs = c("cc","cc"), k=c(40,15)), 
                data=bet[bet$Site == "Bratt",], family=gaussian,
                method="REML", select=TRUE)
tempmodflex <- gam(temp ~ s(Year, k=5) + te(DayofYear, HourEnding, bs = c("cc","cc"),
                                            by=factor(Year)), 
               data=bet[bet$Site == "Bratt",], family=gaussian,
               method="REML", select=TRUE)

## get expand grid object to plot a ggplot contour/heatmap
preddf <- with(bet, expand.grid(HourEnding = seq(min(HourEnding), max(HourEnding), length = 300), 
                                     DayofYear = seq(min(DayofYear), max(DayofYear), length = 300)))
superN <- nrow(preddf)
YearMean <- with(bet, mean(Year))
preddf$Year <- YearMean

comb.pred <- predict(brattmod, newdata = preddf, type = "terms")

whichCols <- grep("Hour", colnames(comb.pred))
comb.pdat <- cbind(preddf, Fitted = comb.pred[, whichCols])

shiftcomb <- attr(comb.pred, "constant")
comb.pdatnorm <- comb.pdat
comb.pdatnorm <- with(comb.pdatnorm, transform(comb.pdatnorm, Fitted = Fitted + shiftcomb))

#exclude things too far from real
toofar <- exclude.too.far(comb.pdatnorm$DayofYear, comb.pdatnorm$HourEnding, 
                          bet$DayofYear, bet$HourEnding, dist=0.1)
comb.pdatnorm$CO2 <- comb.pdatnorm$Fitted
comb.pdatnorm$CO2[toofar] <- NA

co2predicted <- comb.pdatnorm

comboplot <- ggplot(comb.pdatnorm, aes(x = DayofYear, y =HourEnding, z=CO2)) + #, z=Fitted
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=CO2)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent', name=expression(paste(italic('p')*'CO'[2]*" (ppm)"))) +
  #geom_point(data=bet, aes(x=DayofYear, y=HourEnding, z=NULL)) +
  theme(legend.key.width=unit(2,"cm")) +
  ylab("Hour") + xlab("Day of Year") +
  guides(guide=guide_legend(title.vjust = 0.1))
## FIXME: this does nothing to change the justification of pCO2!!

## for temp:
preddf <- with(bet, expand.grid(HourEnding = seq(min(HourEnding), max(HourEnding), length = 300), 
                                DayofYear = seq(min(DayofYear), max(DayofYear), length = 300)))
superN <- nrow(preddf)
YearMean <- with(bet, mean(Year))
preddf$Year <- YearMean

comb.pred <- predict(tempmod, newdata = preddf, type = "terms")

whichCols <- grep("Hour", colnames(comb.pred))
comb.pdat <- cbind(preddf, Fitted = comb.pred[, whichCols])

shiftcomb <- attr(comb.pred, "constant")
comb.pdatnorm <- comb.pdat
comb.pdatnorm <- with(comb.pdatnorm, transform(comb.pdatnorm, Fitted = Fitted + shiftcomb))

#exclude things too far from real
toofar <- exclude.too.far(comb.pdatnorm$DayofYear, comb.pdatnorm$HourEnding, 
                          bet$DayofYear, bet$HourEnding, dist=0.1)
comb.pdatnorm$temp <- comb.pdatnorm$Fitted
comb.pdatnorm$temp[toofar] <- NA


tempplot <- ggplot(comb.pdatnorm, aes(x = DayofYear, y =HourEnding, z=temp)) + #, z=Fitted
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=temp)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent', name=expression(paste("temperature"*'CO'[2]*" (ppm)"))) +
  #geom_point(data=bet, aes(x=DayofYear, y=HourEnding, z=NULL)) +
  theme(legend.key.width=unit(2,"cm")) +
  ylab("Hour") + xlab("Day of Year") +
  guides(guide=guide_legend(title.vjust = 0.1))
## model seems to need an interaction for year and day of year ... 

## for temp by year:
preddf <- with(bet, expand.grid(HourEnding = seq(min(HourEnding), max(HourEnding), length = 300), 
                                DayofYear = seq(min(DayofYear), max(DayofYear), length = 300)))
superN <- nrow(preddf)
YearN <- length(unique(bet$Year))
YearZ <- unique(bet$Year)
combdf <- data.frame(HourEnding = rep(preddf$HourEnding, times=YearN),
                     DayofYear = rep(preddf$DayofYear, times=YearN),
                     Year = rep(YearZ, each=superN))

comb.pred <- predict(tempmodflex, newdata = combdf, type = "terms")

whichCols <- grep("Hour", colnames(comb.pred))
comb.pdat <- cbind(preddf, Fitted = rowSums(comb.pred[, whichCols]))
comb.pdat$Year <- rep(YearZ, each=superN)

shiftcomb <- attr(comb.pred, "constant")
comb.pdatnorm <- comb.pdat
comb.pdatnorm$Fitted <- comb.pdatnorm$Fitted + shiftcomb

#exclude things too far from real
toofar <- exclude.too.far(comb.pdatnorm$DayofYear, comb.pdatnorm$HourEnding, 
                          bet$DayofYear, bet$HourEnding, dist=0.1)
comb.pdatnorm[toofar,] <- NA

tempplotflex <- ggplot(comb.pdatnorm, aes(x = DayofYear, y =HourEnding, z=Fitted, group=Year)) + 
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=Fitted)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent', name=expression(paste("Fill: Temperature ("~degree*"C)"))) +
  #geom_point(data=bet, aes(x=DayofYear, y=HourEnding, z=NULL)) +
  theme(legend.key.width=unit(2,"cm")) +
  facet_wrap("Year") +
  ylab("Hour") + xlab("Day of Year") +
  geom_contour(inherit.aes = FALSE, aes(x=DayofYear, y=HourEnding, z=CO2, colour=..level..))+
  scale_color_viridis(name=expression(paste("Lines: Multiyear mean"~italic('p')*'CO'[2]*" (ppm)")))


  grid.arrange(comboplot, tempplot)

## save plots and models
saveRDS(brattmod, "../data/private/bratt-gam-co2time.rds")
saveRDS(tempmod, "../data/private/bratt-gam-temptime.rds")
saveRDS(tempmodflex, "../data/private/bratt-gam-temptime-byyear.rds")
ggsave("../docs/private/bratt-co2station.pdf", brattplot, scale=1.6)
