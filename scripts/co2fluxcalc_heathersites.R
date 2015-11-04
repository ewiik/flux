## run gasexchange on heather's sites

## load required data
if (!file.exists("data/private/heathergasfluxsupp.rds")) {
  source("scripts/co2flux_heathersites.R")
}
alldat <- readRDS("data/private/heathergasfluxsupp.rds")
alldat <- transform(alldat, sampleDate = as.POSIXct(as.character(sampleDate), format = "%Y-%m-%d"))
colnames(alldat)[which(colnames(alldat) == "sampleDate")] <- "Date"
alldat <- transform(alldat, Year = as.numeric(format(Date, format = "%Y")),
                    Month = as.numeric(format(Date, format = "%m")),
                    Day = as.numeric(format(Date, format = "%d")))

ml <- read.csv("data/maunaloa.csv") 

## Insert pco2atm from Mauna Loa
mlsub <- subset(ml, select = c('year', 'month', 'value'))
names(mlsub) <- c("Year", "Month", "pco2atm")
alldat <- merge(alldat, mlsub, by = c("Year", "Month"))

## load gasExchangeFlex
source("functions/gasExchangeFlex.R")

## run the function, feed in changing atmospheric co2 concs
args(gasExchangeFlex) # couldn't remember...
alldat <- transform(alldat, co2flux = gasExchangeFlex(temp = temperature, cond = conductivity, ph = pH,
                                                      wind = wind, salt = salinity, dic = DIC, 
                                                      alt = Altitude, altnotkpa = TRUE,
                                                      pco2atm = pco2atm))

## a few diagnostic plots...
with(alldat, plot(co2flux ~ salinity))
with(alldat, plot(co2flux ~ conductivity))
with(alldat, plot(co2flux ~ Altitude))
with(alldat, plot(co2flux ~ wind))
with(alldat, plot(co2flux ~ temperature))
with(alldat, plot(co2flux ~ DIC))
with(alldat, plot(co2flux ~ pH))
with(alldat, plot(latitude ~ longitude, type = "n"))
with(alldat, points(latitude ~ longitude, cex = abs(co2flux)/10))
with(alldat, points(latitude ~ longitude, cex = abs(salinity/5)))
# could also do as per symbols(x=dfx$ev1, y=dfx$ev2, circles=dfx$ev3, inches=1/3,
# ann=F, bg="steelblue2", fg=NULL)
