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

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=18, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
esther <- read.delim("~/not-to-backup/Esther-CO2-Hourly.DAT", sep = ",")
bratt <- read.delim("~/not-to-backup/Bratts_Lake-CO2-Hourly.DAT", sep = ",")
trout <- read.delim("~/not-to-backup/East_Trout_Lake-CO2-Hourly.DAT", sep = ",")

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

bet <- merge(bet, sundat) 
bet$TOD <- ifelse(bet$HourEnding-0.5 < bet$DownTime & 
                    bet$HourEnding-0.5 > bet$UpTime, 
                  "Day", "Night")

## plot them and look
brattplot <- ggplot(data=bet[bet$Year>=2009 & bet$Site == "Bratt",], aes(y=Mean, x=DayofYear, col=TOD)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_wrap("Year", ncol=4) +
  scale_color_viridis(name="Time of Day", discrete=TRUE, option = "magma", direction = -1,
                      begin=0, end=0.7)  +
  #scale_color_viridis(name="Hour") + for when plotting by Hour
  theme(legend.key.width=unit(2,"cm")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylab("CO2 (ppm)") + xlab("Day of Year (200 mid-July)")

ggplot(data=bet[bet$Year>=2009 & bet$Site == "Esther",], aes(y=Mean, x=DayofYear, col=TOD)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_wrap("Year", ncol=4) +
  scale_color_viridis(name="Time of Day", discrete=TRUE, option = "magma", direction = -1,
                      begin=0, end=0.7)  +
  #scale_color_viridis(name="Hour") + for when plotting by Hour
  theme(legend.key.width=unit(2,"cm")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  ylab("CO2 (ppm)") + xlab("Day of Year (200 mid-July)")

ggplot(data=bet, aes(y=Mean, x=DayofYear, col=TOD)) +
  papertheme +
  geom_point(alpha=0.3) +
  facet_grid(Year ~ Site) +
  scale_color_viridis(name="Time of Day", discrete=TRUE, option = "magma", direction = -1,
                      begin=0, end=0.7)  +
  #scale_color_viridis(name="Hour") + for when plotting by Hour
  theme(legend.key.width=unit(2,"cm")) +
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
brattmod <- gam(Mean ~ s(Year, k=5) + te(DayofYear, HourEnding, bs = c("cc","cc"), k=c(20,15)), 
                data=bet[bet$Site == "Bratt",], family=scat(),
                method="REML", select=TRUE)
## save plots and models
ggsave("../docs/private/bratt-co2station.pdf", brattplot, scale=1.6)
