## This script creates rds's that contain all relevant data for running regressions against pCO2:
## pdo.rds, soi_nonstand.rds, soi_stand.rds, naoseasonal.rds, co2explained.rds
## 1. download weather pattern data from web: PDO, SOI, NAO
## 2. process routines sampling data tables
## 3. get flow et al data from other sources

## Part 1: all these data are complete
## ==================================================================================================

## download PDO url and create data frame structure with numerics
if (!file.exists("data/pdo.txt")) { 
  download.file("http://research.jisao.washington.edu/pdo/PDO.latest", "data/pdo.txt")
  }

file <- read.fwf(file = "data/pdo.txt", skip = 30, n = 118, 
                   widths = c(8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7))

colnames(file) <- unlist(file[1,])
pdo <- file[-c(1, 2),] # all as factors and straight to numeric gives gibberish
pdo <- lapply(pdo, as.character)
pdo <- lapply(pdo, as.numeric) #though this makes NA out of the years that have asterisks.
#   these asterisks actually mean "Derived from OI.v2 SST fields" (erm...)
pdo <- do.call(cbind, pdo)
pdo[pdo == 99.9] <- NA
pdo <- as.data.frame(pdo)
names(pdo) <- c("year", "jan", "feb", "mar", "apr","may", "jun", "jul", "aug", "sep", "oct", "nov", 
                "dec") # unlist() earlier maintained spaces in the names...
# could've had month.abb here to c all months in a year
allna <- which(is.na(pdo$year))
pdo$year[allna] <- c(seq(from = 2002, length.out = length(allna), by = 1))

## create yearly means since I believe we want annual reso
pdo <- transform(pdo, mean = rowMeans(pdo[,-1], na.rm = TRUE))

## create rds for later
saveRDS(pdo, "data/pdo.rds")

## download SOI url and create data frames with numeric entries
## FIXME: do we want standardised or unstandardised data?
if (!file.exists("data/soi.txt")) { 
  download.file("http://www.cpc.noaa.gov/data/indices/soi", "data/soi.txt")
  }

all_data <- read.fwf(file = "data/soi.txt", skip = 3, widths = rep(6, 13))

## create separate df for the unstandardised data
table1 <- all_data[c(1:66),] 
colnames(table1) <- unlist(table1[1,])
soinonst <- table1[-c(1),] # all as factors and straight to numeric gives gibberish
soinonst <- lapply(soinonst, as.character)
soinonst <- lapply(soinonst, as.numeric)
soinonst <- do.call(cbind, soinonst)
soinonst[soinonst == 99.9] <- NA
soinonst <- as.data.frame(soinonst)

## create rds for later
saveRDS(soinonst, "data/soi_nonstand.rds")

## create separate df for the standardised data
table2 <- all_data[c(75:140),] 
colnames(table2) <- unlist(table2[1,])
soist <- table2[-c(1),] # all as factors and straight to numeric gives gibberish
soist <- lapply(soist, as.character)
soist <- lapply(soist, as.numeric)
soist <- do.call(cbind, soist)
soist[soist == 99.9] <- NA
soist <- as.data.frame(soist)

## create rds for later
saveRDS(soist, "data/soi_stand.rds")


## download NAO url and make data frame with numeric entries
if (!file.exists("data/nao.txt")) {
  download.file("https://climatedataguide.ucar.edu/sites/default/files/nao_station_seasonal.txt",
                "data/nao.txt")
}

naos <- readLines("data/nao.txt")
naos <- naos[-c(1:2)]
naos <- lapply(naos, gsub, pattern = "    ", replacement = " ")
naos <- lapply(naos, gsub, pattern = "   ", replacement = " ")
naos <- lapply(naos, gsub, pattern = "  ", replacement = " ")

naosplit <- lapply(naos, strsplit, " ")
naovec <- unlist(naosplit)
naomat <- matrix(naovec, ncol=13, byrow=TRUE)
naoframe <- as.data.frame(naomat)
colnames(naoframe) <- c("year", "djf", "jfm", "fma", "mam", "amj", "mjj", "jja", "jas", "aso", 
                      "son", "ond", "ndj")
naoframe <- lapply(naoframe, as.character)
naoframe <- lapply(naoframe, as.numeric)
naoframe <- do.call(cbind, naoframe)
naoframe <- as.data.frame(naoframe)
naoframe[naoframe <= -999] <- NA

## create rds for later
saveRDS(naoframe, "data/naoseasonal.rds")

## Part 2: 1039 rows, 415 rows of one or more NAs 
## ==================================================================================================

## read in supporting monitoring data tables grabbed from database
## I created a Date2 column in OpenOffice to convert the current display of month as character to month
##    as numeric (former is how it came from database into excel)

routines <- read.csv("data/private/qpco2supportdata.csv")
routines <- transform(routines, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))
routines <- routines[with(routines, order(LAKE, Date)),]

produc <- read.csv("data/private/Production.csv")
produc <- transform(produc, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))

chl <- read.csv("data/private/Chl.csv")
chl <- transform(chl, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))

incub <- read.csv("data/private/qprodsupportdata.csv")
incub <- transform(incub, Date = as.POSIXct(as.character(Date), format = "%d-%m-%Y"))

oxtemp <- read.csv("data/private/qprofilesoxtemp.csv")
oxtemp <- rbind(oxtemp, read.csv("data/private/qprofilesextrarecordsoxtemp.csv"))
oxtemp <- transform(oxtemp, Date = as.POSIXct(as.character(Date), format = "%d-%m-%Y"))


secchi <- read.csv("data/private/qsecchietal.csv") 
secchi <- transform(secchi, Date = as.POSIXct(as.character(Date), format = "%Y-%m-%d"))

lakes <- read.csv("data/private/Lakes.csv") 
colnames(lakes)[which(colnames(lakes) == "Abbreviation")] <- "LAKE"

## remove extra date column (originally retained in case wanna check that the date format
##    conversion worked in OpenOffice)
routines <- routines[,-which(colnames(routines)=="Date2")] 
produc <- produc[,-which(colnames(produc)=="Date2")] 
chl <- chl[,-which(colnames(chl)=="Date2")] 

## !!!!!!! FIXME: In database: lake names other than WW/WC *at least* in 2012 have been entered
##    with a trailing space which means that merges won't work!
chl <- transform(chl, LAKE = gsub(" ", "", LAKE))


## take out all chl data that isn't from an integrated sample, as that's what we're working with
##    (though worth thinking that surface could be played with since they're quite different!!)
## FIXME: All 63micron samples are from 1997 - what is this sample, and should it be included?
chlsub <- subset(chl, TreatmentNewLabel == "Integrated")

## split produc and chl for getting means for replicated production estimates
prodsplit <- with(produc, split(produc, list(LAKE, Date), drop = TRUE))
chlsplit <- with(chlsub, split(chlsub, list(LAKE, Date), drop = TRUE))

## create wrapper that does colmeans for the columns we want and spits out df we can merge
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(LAKE = df['LAKE'][1,], Date = df['Date'][1,]), t(means))
}

## choose columns we want means for
## FIXME: chl columns: since we take "chl a" rather than total chl there are NAs; should we try 
##    total chl too??

chlmeans <- lapply(chlsplit, mycolMeans, cols = c("Chl_a_ug_L")) # , "Total_chl"
chlmeans <- do.call(rbind, chlmeans)
rownames(chlmeans) <- NULL
chlmeans <- chlmeans[with(chlmeans, order(LAKE, Date)),]

prodmeans <- lapply(prodsplit, mycolMeans, cols = c("LIGHT_O2_ppm", "DARK_O2_ppm", "LIGHT_Pois_O2_ppm",
                                                    "DARK_Pois_O2_ppm", "NetOxy_ppm", "RespOxy_ppm",
                                                    "Net_mgC_m3_h", "Resp_mgC_m3_H", 
                                                    "LightAzidemgCm3h", "DarkAzidemgCm3h")) 
prodmeans <- do.call(rbind, prodmeans)
rownames(prodmeans) <- NULL
prodmeans <- prodmeans[with(prodmeans, order(LAKE, Date)),]

## merge prodmeans with the other prod support data in another table
## using AlainsCode.R in private/ which nicole sent me, to estimate P's and R.
## FIXME: why are the poisoned ones not used for anything? Nicole says not used but she doesn't
##    know why either; C columns in prodmeans were calculated based on a metabolic ratio 
proddf <- merge(prodmeans, incub, by = c("LAKE", "Date"))
proddf <- transform(proddf, photO2ppm = LIGHT_O2_ppm - ProdAssayO2Start_ppm)
proddf <- transform(proddf, respO2ppm = DARK_O2_ppm - ProdAssayO2Start_ppm)
proddf <- transform(proddf, NPP_h = photO2ppm/Hours_incubation)
proddf <- transform(proddf, R_h = respO2ppm/Hours_incubation)
proddf <- transform(proddf, GPP_h = NPP_h - R_h)

## leave proddf for later if need to revisit calcs, select what we want
prodsub <- subset(proddf, select = c(LAKE, Date, GPP_h, NPP_h, R_h))

### nicole noticed an outlier and there it is: WW 2010-05-03 (GPP >15 cf mostly <1)
## FIXME: changing this to NA for now but may want to keep it and actually deal with all outliers
##    in one sitting?
## FIXME: is the outlier just a typo? could follow up
makena <- which(prodsub$GPP_h > 15)
prodsub[makena, c('GPP_h', 'NPP_h', 'R_h')] <- NA

## finlay et al lakewide estimate method doesn't make sense so doing this as alain and nicole have done
##    (finlay et al imply that secchi divided by lake volume, but this gives a strange proportion rather
##    than an upscaling to whole-lake NPP)
prodsub <- merge(prodsub, lakes[,c('LAKE', 'LakeArea_km2')])
prodsub <- merge(prodsub, secchi[c('LAKE', 'Date', 'Secchi_m')])
prodsub <- transform(prodsub, lakeGPP = GPP_h*(Secchi_m * (LakeArea_km2 * 100000)))
prodsub <- transform(prodsub, lakeNPP = NPP_h*(Secchi_m* (LakeArea_km2 * 100000)))
prodsub <- transform(prodsub, lakeR = R_h*(Secchi_m* (LakeArea_km2 * 100000)))

## merge all dfs by LAKE and Date
co2explained <- merge(chlmeans, prodsub, all = TRUE) ## FIXME: chl NAs are now NaN.. problem for later?
co2explained <- merge(co2explained, routines[,-which(names(routines) %in% 'RunNo')], all = TRUE) 
# don't need this column which is cryptic anyway
co2expl <- merge(co2explained, oxtemp[,c('Date', 'LAKE', 'Temperature_deg_C','Oxygen_ppm')], 
                 by = c('Date', 'LAKE'), all = TRUE)
## FIXME: minor discrepancies in nrow over merge steps... a few rows.. could check at some point

## airtemp, from env canada and weathermanipulations.R
if (!file.exists("data/temperaturedata.rds")) {
  source("scripts/weathermanipulations.R")
}
airtemp <- readRDS("data/temperaturedata.rds")
## since we decided to use the regina station for all measurements, I will also here subset
##    to regina
airtemp <- subset(airtemp, Superstation == "regina", select = c(Year, Month, Temperature))
airtemp <- airtemp[order(airtemp$Year, airtemp$Month),]
rownames(airtemp) <- NULL #set new order for rows, otherwise in mixed order based on above command

## create yearmeans for temperature, too
airsplit <- with(airtemp, split(airtemp, list(Year)))
airmeans <- lapply(airsplit, colMeans, na.rm = TRUE)
airmeans <- do.call(rbind, airmeans)
colnames(airmeans)[which(colnames(airmeans) == "Temperature")] <- "AirTempAnnual"
airmeans <- airmeans[,-(which(colnames(airmeans) == "Month"))]
airmeans <- as.data.frame(airmeans)
airmeans$Year <- as.numeric(airmeans$Year)
airtemp <- merge(airtemp, airmeans)
colnames(airtemp)[which(colnames(airtemp) == "Temperature")] <- "AirTempMonthly"
colnames(airtemp)[which(colnames(airtemp) == "Year")] <- "YEAR"

## merge co2expl with monthly air temperature data....
co2expl <- transform(co2expl, Month = as.numeric(format(Date, format = "%m")))
co2expl <- merge(co2expl, airtemp)

## relhum, from env canada and weathermanipulations.R
if (!file.exists("data/relhumdata.rds")) {
  source("scripts/weathermanipulations.R")
}
relhum <- readRDS("data/relhumdata.rds")
## since we decided to use the regina station for all measurements, I will also here subset
##    to regina
relhum <- subset(relhum, Superstation == "regina", select = c(Year, Month, RelHum))
relhum <- relhum[order(relhum$Year, relhum$Month),]
names(relhum)[which(names(relhum) == "Year")] <- "YEAR"
rownames(relhum) <- NULL #set new order for rows, otherwise in mixed order based on above command

## merge co2expl with monthly relative humidity data
co2expl <- merge(co2expl, relhum)

## save output for later
saveRDS(co2expl, "data/private/co2explained.rds")

## bonus: let's see how many NAs in our data. (don't select Si and other optionals)
dataloss <- subset(co2expl, select = c("LAKE", "Date", "pH_surface", "Oxygen_ppm", "TIC_mg_L", 
                                       "Temperature_deg_C", "SRP_ug_L", "TDN_ug_L", "NO3_ug_L", 
                                       "NH4_ug_L", "DOC_mg_L", "Chl_a_ug_L", "GPP_h"))

nanumbers <- rowSums(is.na(dataloss))
dataloss <- dataloss[rowSums(is.na(dataloss)) > 0,]
nrow(dataloss)

## Part 3: ice-out, inflow, evaporation
## ==========================================================================================

## kerri emailed me this (in PhDPapers/Data..) pCO2_vars_yearlyaverageswithclimate.csv, but
##    it has only a few years, so we need to supplement it.
## FIXME: need to know which flow site for which lake
## other files downloaded at this point are from Rich. Go to 2010/2011
## kerri used the "old" way of calculating annual flows, but nicole and I checking what data
##    we can get to bring us up to 2015. probably will be the "new way"
## "new way": monthly data set from rich (full csv in !git/fromrich, sensible csv in /private..)
##    filled in a few missing annual totals in openoffice before saving as csv
## "old way": annual data set from rich, full csv in !git/fromrich, sensible csv in /private..)
##  Interestingly: why does finlay et al 2009 say that inflow for crooked was regressed, since it 
##    exists in the data set that rich sent that should correspond to what kerri did?

## ice-out
iceout1 <- read.csv("data/private/IceOut.csv") # pending data to 2015
icedataloss <- iceout1[is.na(iceout1$ICEOUTDOY),]
nrow(iceout1) # 120
nrow(icedataloss) # 45

## inflow flows
inflow1 <- read.csv("data/private/Monthly_Inflow_Estimates2.csv") # pending data to 2015
inflowold <- read.csv("data/private/Lake_Inflows.csv") # 1994-2009

othervars <- read.csv("data/private/pCO2_vars_yearlyaverageswithclimate.csv")
## this doesn't however indicate which vars interpolated and which real, but has the separate inflow
##    sites and inflow data
evap <- subset(othervars, select = "evap" %in% colnames(othervars))
