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

## Part 2: 939 rows, 375 rows of one or more NAs with variable selection FIXMEs 
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

gasflux <- readRDS("data/private/gasFlux.rds") # from pressuremanipulations.R

## remove extra date column (originally retained in case wanna check that the date format
##    conversion worked in OpenOffice)
routines <- routines[,-which(colnames(routines)=="Date2")] # wanted to avoid using numbers for columns
produc <- produc[,-which(colnames(produc)=="Date2")] 
chl <- chl[,-which(colnames(chl)=="Date2")] 

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
## FIXME: confirm with kerri that we want netoxy, and which respiration measure we want; NP vs GPP?
##    apparently this may be a macro that alain calculated in systat then R (see code that Nicole 
##    commented and sent: AlainsCode.R in private/)
## FIXME: is the variable O2 from the water column or the netoxy thing???
## FIXME: chl columns: if we take "chl a" rather than total chl there are NAs --> implies different
##    method used for total chl and individual chlorophylls...!?
##    so need to know if we are taking total or chl a ... and if we can replace missing data

chlmeans <- lapply(chlsplit, mycolMeans, cols = c("Chl_a_ug_L", "Total_chl")) 
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
proddf <- merge(prodmeans, incub, by = c("LAKE", "Date"))



## merge all dfs by LAKE and Date
co2explained <- merge(chlmeans, prodmeans) ## FIXME: NAs are now NaN.. problem for later?
co2explained <- merge(co2explained, routines)
names(gasflux)[which(names(gasflux) == 'Lake')] <- "LAKE"
co2expl <- merge(co2explained, gasflux[,c('Date', 'LAKE', 'Temperature')], by = c('Date', 'LAKE'))
# this creates shorter data frame cause co2explained had R, E, M, WC in it and gasflux not. though 
#   FIXME: diff of 100 rows indicated by nrow() of the dfs vs 97 rows in co2explained when subsetting
#   nrow(subset(co2explained, LAKE %in% c("R", "E", "WC", "M")))

## save output for later
saveRDS(co2expl, "data/private/co2explained.rds")

## bonus: let's see how many NAs in our data. though see FIXME about variable selection
dataloss <- subset(co2expl, select = c("LAKE", "Date", "pH_surface", "TIC_mg_L", "Temperature", 
                                       "SRP_ug_L", "TDN_ug_L", "NO3_ug_L", "NH4_ug_L", "DOC_mg_L", 
                                       "Chl_a_ug_L", "NetOxy_ppm", "Resp_mgC_m3_H"))

nanumbers <- rowSums(is.na(dataloss))
dataloss <- dataloss[rowSums(is.na(dataloss)) > 0,]
 

## Part 3:
## ==========================================================================================

## kerri emailed me this (in PhDPapers/Data..) pCO2_vars_yearlyaverageswithclimate.csv, but
##    it has only a few years, so we need to supplement it.
## flow data summaries can be sourced from https://www.wsask.ca/Lakes-and-Rivers/
##                               Stream-Flows-and-Lake-Levels/QuAppelle-River-Watershed/05JK007/
## but FIXME: need to know which site for which lake

## ice-out, some flow data, from Rich (.../fromRich) and copied into git../private/. 
##    but generally only to 2009-2011. very patchy data

iceout1 <- read.csv("data/private/IceOut.csv") # pending data all the way to 2015
icedataloss <- iceout1[rowSums(is.na(iceout1)) > 0,]
nrow(iceout1)
nrow(icedataloss)

inflow1 <- read.csv("data/private/Monthly_Inflow_Estimates2.csv") # csv not currently in good format
#   but just wanted to check how complete - data are complete.

othervars <- read.csv("data/private/pCO2_vars_yearlyaverageswithclimate.csv")
## this doesn't however indicate which vars interpolated and which real, but has the separate inflow
##    sites and inflow data


