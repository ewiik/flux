## run gasexchange on heather's sites
## chl data is from the database too, heather emailed it on Jan30-2016
## FIXME: Which Little Manitou is in the actual chl data?

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


## Insert pco2atm from Mauna Loa
mlsub <- subset(ml, select = c('Year', 'Month', 'pCO2'))
names(mlsub) <- c("Year", "Month", "pco2atm")
alldat <- merge(alldat, mlsub, by = c("Year", "Month"))

## load gasExchangeFlex
source("../functions/gasExchangeFlex.R")

## run the function, feed in changing atmospheric co2 concs
args(gasExchangeFlex) # couldn't remember...
co2z <- with(alldat, gasExchangeFlex(temp = temperature, cond = conductivity, ph = pH,
                          wind = wind, salt = salinity, dic = DIC, 
                          alt = Altitude, altnotkpa = TRUE,
                          pco2atm = pco2atm))

alldatz <- cbind(alldat, co2z)

## a few diagnostic plots...
pdf("data/private/heathersites.pdf")
with(alldat, plot(co2flux ~ salinity))
with(alldat, plot(co2flux ~ conductivity))
with(alldat, plot(co2flux ~ Altitude))
with(alldat, plot(co2flux ~ wind))
with(alldat, plot(co2flux ~ temperature))
with(alldat, plot(co2flux ~ DIC))
with(alldat, plot(co2flux ~ pH))
with(alldat, plot(latitude ~ longitude, type = "p", pch = 25, bg = "black", cex = 0.4))
with(alldat, points(latitude ~ longitude, cex = abs(co2flux)/10, 
                    col = ifelse(co2flux < 0, "black", "purple")))
with(alldat, points(latitude ~ longitude, cex = abs(salinity/5), col = "green"))
legend("bottomright", cex = 0.8, legend = c("lake location", "salinity", "CO2 influx"), 
       pch = c(25,1,1), pt.bg = c("black", "white", "white"), col = c("black", "green", "black"))
# could also do as per symbols(x=dfx$ev1, y=dfx$ev2, circles=dfx$ev3, inches=1/3,
# ann=F, bg="steelblue2", fg=NULL)
dev.off()

## subset by depth <= 2.5 for kerri
puddles <- subset(alldat, select = c("Date", "lakeName", "latitude", "longitude", "sampleDepth",
                                     "Altitude", "co2flux"), subset = sampleDepth <= 2.5)
write.csv(puddles, "data/private/co2fluxheathersubset.csv")


## =====================================================================================
## look at chlorophyll data in relation with gas flux.
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

alldatchl <- merge(alldatz, chlormeans, all.x = TRUE, by.x = 'lakeName', by.y = 'LAKE')

write.csv(alldatchl, "../data/private/co2fluxheathersites.csv")
