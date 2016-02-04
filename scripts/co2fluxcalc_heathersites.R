## run gasexchange on heather's sites
## chl data is from the database too, heather emailed it on Jan30-2016
## Little Manitou in the chl data set is the north one, site id 66
## FIXME: complete units: secchi & sample depth (m), wind (m/s), conductivity (uS/cm), .. 
## FIXME: incorporate these into final colnames

## load required data
if (!file.exists("../data/private/heathergasfluxsupp.rds")) {
  source("co2flux_heathersites.R")
}
alldat <- readRDS("../data/private/heathergasfluxsupp.rds")
alldat <- transform(alldat, sampleDate = as.POSIXct(as.character(sampleDate), format = "%Y-%m-%d"))
colnames(alldat)[which(colnames(alldat) == "sampleDate")] <- "Date"
alldat <- transform(alldat, Year = as.numeric(format(Date, format = "%Y")),
                    Month = as.numeric(format(Date, format = "%m")),
                    Day = as.numeric(format(Date, format = "%d")))

ml <- read.csv("../data/maunaloa.csv") 
if (!file.exists("../data/private/chlorophyllheathersites.csv")) {
  stop("get chlorophyll data from heather")
} else {chlor <- read.csv("../data/private/chlorophyllheathersites.csv")}

## load required function
source("../functions/gasExchangeFlex.R")

## Insert pco2atm from Mauna Loa
mlsub <- subset(ml, select = c('Year', 'Month', 'pCO2'))
names(mlsub) <- c("Year", "Month", "pco2atm")
alldat <- merge(alldat, mlsub, by = c("Year", "Month"))

## convert DIC from mg/L to umol
alldat$DICumol <- alldat$DIC / 0.012

## add chlorophyll data for comparison with gas flux.
chlor <- subset(chlor, chlID == "A")
chlormeans <- with(chlor, split(chlor, list(sampleID)))

mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(LAKE = df['sampleID'][1,], t(means)))
}

chlormeans <- lapply(chlormeans, mycolMeans, 'chlvalue')
chlormeans <- do.call(rbind, chlormeans)
row.names(chlormeans) <- NULL

alldat <- merge(alldat, chlormeans, all.x = TRUE, by.x = 'lakeName', by.y = 'LAKE')

## run the function, feed in changing atmospheric co2 concs
args(gasExchangeFlex) # couldn't remember...
co2z <- with(alldat, gasExchangeFlex(temp = temperature, cond = conductivity, ph = pH,
                          wind = wind, salt = salinity, dic = DICumol, 
                          alt = Altitude, altnotkpa = TRUE,
                          pco2atm = pco2atm))

alldatz <- cbind(alldat, co2z)
names(alldatz)[which(names(alldatz) %in% c("fluxenh", "pco2"))] <-
  c("CO2fluxmmolm2d", "pCO2uatm")

## save object for later
write.csv(alldatz, "../data/private/co2fluxheathersites.csv", row.names = FALSE)

## a few diagnostic plots...
pdf("data/private/heathersites.pdf")
with(alldatz, plot(CO2fluxmmolm2d ~ salinity))
with(alldatz, plot(CO2fluxmmolm2d ~ conductivity))
with(alldatz, plot(CO2fluxmmolm2d ~ Altitude))
with(alldatz, plot(CO2fluxmmolm2d ~ wind))
with(alldatz, plot(CO2fluxmmolm2d ~ temperature))
with(alldatz, plot(CO2fluxmmolm2d ~ DIC))
with(alldatz, plot(CO2fluxmmolm2d ~ pH))
with(alldatz, plot(latitude ~ longitude, type = "p", pch = 25, bg = "black", cex = 0.4))
with(alldatz, points(latitude ~ longitude, cex = abs(CO2fluxmmolm2d)/10, 
                    col = ifelse(CO2fluxmmolm2d < 0, "black", "purple")))
with(alldatz, points(latitude ~ longitude, cex = abs(salinity/5), col = "green"))
legend("bottomright", cex = 0.8, legend = c("lake location", "salinity", "CO2 influx"), 
       pch = c(25,1,1), pt.bg = c("black", "white", "white"), col = c("black", "green", "black"))
# could also do as per symbols(x=dfx$ev1, y=dfx$ev2, circles=dfx$ev3, inches=1/3,
# ann=F, bg="steelblue2", fg=NULL)
dev.off()

with(alldatz, plot(pCO2uatm ~ salinity))
with(alldatz, plot(pCO2uatm ~ log(conductivity)))
with(alldatz, plot(pCO2uatm ~ Altitude))
with(alldatz, plot(pCO2uatm ~ wind))
with(alldatz, plot(pCO2uatm ~ temperature))
with(alldatz, plot(pCO2uatm ~ DIC))
with(alldatz, plot(pCO2uatm ~ pH))
abline(0,0)
with(alldatz, plot(pCO2uatm ~ chlvalue, xlim = c(0,200))))


## =============================================================================================
## subset by depth <= 2.5 for kerri
## =============================================================================================
puddles <- subset(alldatz, select = c("Date", "lakeName", "latitude", "longitude", "sampleDepth",
                                     "Altitude", "co2flux"), subset = sampleDepth <= 2.5)
write.csv(puddles, "data/private/co2fluxheathersubset.csv")

