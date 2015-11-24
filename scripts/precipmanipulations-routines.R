## use gasexchange_precipdata.R to source regina info for routines analysis

## source necessary functions:
source("functions/gasexchange_precipdata.R")

## my station data frame: for regina gilmour as per heather's searches
stations <- data.frame(StationID = c(3007), start = c(1994), end = c(2014))

precip <- getData(stations, folder = "data/precipdata/routines")
saveRDS(precip, "data/precip.rds") # save for future
