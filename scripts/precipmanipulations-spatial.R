## use gasexchange_precipdata.R to source weather station info for spatial analysis:
##    heather's sites

## source necessary functions:
source("functions/gasexchange_precipdata.R")

## heather's station data frame
stations <- data.frame(StationID = c(), # insert list of station IDs here
                       start = c(), # insert first year of required data for each station                                                                                data start for each station]), 
                       end = c()) # insert last year of required data for each station

precipSK <- getData(stations, folder = "data/precipdata/spatial")

## ANY STATIONS IN ALBERTA??
