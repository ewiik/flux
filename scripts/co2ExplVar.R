## This script creates an rds that contains all relevant data for running regressions against pCO2
## 1. download weather pattern data from web: PDO, SOI, NAO
## 2. process routines sampling data tables
## 3. get flow et al data from other sources

## Part 1:
## ==================================================================================================

## urls, and filenames to use to save the data

dload <- readLines("http://research.jisao.washington.edu/pdo/PDO.latest", n = 148)
dload2 <- dload[-c(1:30, 32)]
dload2[117] <- paste(dload2[117], " NA NA NA NA")
dtable <- dload2

dtable2 <- lapply(dtable, gsub, pattern = "      ", replacement = " ")
dtable2 <- lapply(dtable2, gsub, pattern = "     ", replacement = " ")
dtable2 <- lapply(dtable2, gsub, pattern = "    ", replacement = " ")
dtable2 <- lapply(dtable2, gsub, pattern = "   ", replacement = " ")
dtable2 <- lapply(dtable2, gsub, pattern = "  ", replacement = " ")

dsplit <- lapply(dtable2, strsplit, " ")
dvec <- unlist(dsplit)
dmat <- matrix(dvec, ncol=13, byrow=TRUE)
dframe <- as.data.frame(dmat)

write.csv(dframe, "data/pdo.csv")

## OR, as Jon found read.fwf:
if (file.exists("data/pdo.txt")) { 
  file <- read.fwf(file = "data/pdo.txt", skip = 30, n = 118, 
                   widths = c(8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7))
} else {
  download.file("http://research.jisao.washington.edu/pdo/PDO.latest", "data/pdo.txt")
  file <- read.fwf(file = "data/pdo.txt", skip = 30, n = 118, 
                   widths = c(8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7))
}

colnames(file) <- unlist(file[1,])
pdo <- file[-c(1, 2),]

if (file.exists("data/soi.txt")) { 
  all_data <- read.fwf(file = "data/soi.txt", skip = 3, widths = rep(6, 13))
} else {
  download.file("http://www.cpc.noaa.gov/data/indices/soi", "data/soi.txt")
  all_data <- read.fwf(file = "data/soi.txt", skip = 3, widths = rep(6, 13))
}
table1 <- all_data[c(1:66),] # unstandardised
colnames(table1) <- unlist(table1[1,])
soinonst <- table1[-c(1),]

table2 <- all_data[c(75:140),] # standardised
colnames(table2) <- unlist(table2[1,])
soist <- table2[-c(1),]
## FIXME: do we want standardised or unstandardised?
## FIXME: check that we can read all data as numeric etc.

## Part 2:
## ==================================================================================================

## read in supporting monitoring data tables grabbed from database
## I created a Date2 column in OpenOffice to convert the current display of month as character to month
##    as numeric (former is how it came from database into excel cause excel sucks)

routines <- read.csv("data/private/qpco2supportdata.csv")
routines <- transform(routines, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))
routines <- routines[with(routines, order(LAKE, Date)),]

produc <- read.csv("data/private/Production.csv")
produc <- transform(produc, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))

chl <- read.csv("data/private/Chl.csv")
chl <- transform(chl, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))

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
## FIXME: confirm with kerri that we want netoxy, and which respiration measure we want
## FIXME: chl columns: if we take "total chl", there are NAs where method only measured chla ...
##    so need to know if we are taking total or chl a ... and if we can replace missing data

prodmeans <- lapply(prodsplit, mycolMeans, cols = c("NetOxy_ppm", "Resp_mgC_m3_H")) 
prodmeans <- do.call(rbind, prodmeans)
rownames(prodmeans) <- NULL
prodmeans <- prodmeans[with(prodmeans, order(LAKE, Date)),]

chlmeans <- lapply(chlsplit, mycolMeans, cols = c("CHLC", "Total_chl")) 
chlmeans <- do.call(rbind, chlmeans)
rownames(chlmeans) <- NULL
chlmeans <- chlmeans[with(chlmeans, order(LAKE, Date)),]

## merge all dfs by LAKE and Date
co2explained <- merge(chlmeans, prodmeans) ## FIXME: NAs are now NaN.. problem for later?
co2explained <- merge(co2explained, routines)

## save output for later
saveRDS(co2explained, "data/private/co2explained.rds")

## Part 3:
## ==========================================================================================

## kerri emailed me this (in PhDPapers/Data..) pCO2_vars_yearlyaverageswithclimate.csv, but
##    it has only a few years, so we need to supplement it.
## flow data can be sourced from https://www.wsask.ca/Lakes-and-Rivers/
##                               Stream-Flows-and-Lake-Levels/QuAppelle-River-Watershed/05JK007/
## but FIXME: need to know which site for which lake

## ice-out, some flow data, from Rich (.../fromRich). but generally only to 2009-2011. very patchy data
