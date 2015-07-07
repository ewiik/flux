## we have 1080 dfs in total; in met, they do not have the first 9 lines
## grab met from gasexchange_pressuredata.R

library("reshape")

if (file.exists("dmet.rds")) {
  dmet <- readRDS("dmet.rds")  
} else {
  dmet <- do.call("rbind", met)
  saveRDS(dmet, "dmet.rds")
}


## get station details
stations <- read.csv("data/weatherstations.csv")

## Combine stations
dmet$superstation[dmet$StationID == 3002 | dmet$StationID == 51441] <- "regina"
dmet$superstation[dmet$StationID == 3062 | dmet$StationID == 44203 | dmet$StationID == 48977] <- "yorkton"
dmet$superstation[dmet$StationID == 3318] <- "outlook"
dmet$superstation[dmet$StationID == 2926 | dmet$StationID == 2925] <- "indhead"
unique(dmet$superstation) # yes all appear, no NAs

names(dmet)[2] <- "datetime" # subset can't seem to handle colnames with a slash...

take <- c("StationID","datetime","Year","Month","Day","Time","superstation","Stn Press (kPa)")
dmet <- dmet[take]
names(dmet)[names(dmet) == "Stn Press (kPa)"] <- "Pressure"

## Coerce datetime to the correct internal type
dmet <- transform(dmet, 
                  datetime = as.POSIXct(as.character(datetime), format = "%Y-%m-%d %H:%M"),
                  superstation = as.factor(superstation))

## Delete a priori-known missing data - find station dates at
  ## http://climate.weather.gc.ca/advanceSearch/searchHistoricData_e.html

dmeh <- dmet
dmeh$datetimepos <-as.POSIXlt(dmeh$datetime) 
leaptest <- function(year) {
  if (year %% 4 != 0) { 
    print("not leap") # remaining contenders move on
  }
  else if (year %% 100 == 0 & year %% 400 != 0) { # exclude those that are divisible by both
    print("not leap")
  }
  else {print("is leap")}
} # --> http://disc.gsfc.nasa.gov/julian_calendar.shtml

dmeh$Day[dmeh$datetimepos$year == 113 & dmeh$datetimepos$yday >= 282 & dmeh$StationID == 3002] # year since 1900
maloutflowchemdates[ (as.POSIXlt(maloutflowchemdates)$mon) >= 05 & (as.POSIXlt(maloutflowchemdates)$mon) <= 09], col = "red")

## Remove NA data
miss <- with(dmet, is.na(Pressure))
dmet <- droplevels(dmet[!miss, ])

## Split apart
spldmet <- with(dmet, split(dmet, list(superstation, Year, Month)))

### Wrapper function
meanPressure <- function(df) {
  meanP <- if (NROW(df) == 0) {
    NA
  } else {
    mean(df[["Pressure"]], na.rm = TRUE)
  }
  with(df, data.frame(Superstation = superstation[1],
                      StationID = StationID[1],
                      Year = Year[1],
                      Month = Month[1],
                      Pressure = meanP))
}

pressure <- na.omit(do.call("rbind", lapply(spldmet, meanPressure)))
rownames(pressure) <- NULL


