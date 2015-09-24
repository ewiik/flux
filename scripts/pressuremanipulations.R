## outputting co2 exchange with all data combined as one table
## get necessary data from database queries

if (file.exists("data/pressuredata.rds")) {
  pressuredata <- readRDS("data/pressuredata.rds")  
} else {
  print("run weathermanipulations.R")
  }

if (file.exists("data/windsdata.rds")) {
  windsdata <- readRDS("data/windsdata.rds")  
} else {
  print("run weathermanipulations.R")
}

if (file.exists("data/private/fluxquery1.csv")) { # DIC, pH, wind
  surfd <- read.csv("data/private/fluxquery1.csv")
} else {
  print("get fluxquery1.csv from dropbox")
}

if (file.exists("data/private/fluxquery2.csv")) { # elevation and lake 
  # abbreviation data
  elevd <- read.csv("data/private/fluxquery2.csv")
} else {
  print("get fluxquery2.csv from dropbox")
}

if (file.exists("data/private/fluxquery3.csv")) { # T, cond, salinity
  tcsd <- read.csv("data/private/fluxquery3.csv")
} else {
  print("get fluxquery3.csv from dropbox")
}

if (file.exists("data/superstationz.csv")) { # links superstations with lakes
  superstations <- read.csv("data/superstationz.csv")
} else {
  print("get superstationz.csv from emma")
}

## are we really removing all NA data? THERE'S LOTS! complete.cases, na.omit
## NO - Leave NA's alone so we can see where they are at this stage
remNA <- FALSE
##
surfset <- surfd
# 
tcsset <- tcsd
## Deal with NAs?
if (remNA) {
  tcsset <- droplevels(na.omit(tcsset))
  surfset <- droplevels(na.omit(surfset))
}

## create date objects for later?
tcsset <- transform(tcsset, Date = as.POSIXct(as.character(Date), format = "%d-%m-%Y"))
surfset <- transform(surfset, Date = as.POSIXct(as.character(Date), format = "%d-%m-%Y"))

## merge the tcsset and surfset data sets on
joined <- merge(tcsset, surfset)

## Simplify some names
names(joined) <- c("Lake", "Date", "Temperature", "Conductivity", "Salinity",
                   "Depth", "pH", "Wind", "TIC", "DOY")

## Merge in elevd height data
joined <- merge(joined, subset(elevd, select = c("Abbreviation", "Elevation_m")),
                by.x = "Lake", by.y = "Abbreviation", sort = FALSE)
names(joined)[names(joined) == "Elevation_m"] <- "Elevation"

## drop all lakes we're not interested in 
drop <- TRUE
if(drop) {
  joined <- subset(joined, Lake %in% c("B", "C", "WW", "D", "K", "L", "P"))
}

## choice of using or not using all superstations
regina <- TRUE

if(regina) {
  joined$Superstation <- rep("regina")
} else {
  joined <- merge(joined, superstations, sort = FALSE)
  }

## merge joined with pressure and wind data
## add Year, Month and Day fields to aid the match
## FIXME: Still need to decide if and what to do with missing pressure data.
joined <- transform(joined, Year = as.numeric(format(Date, format = "%Y")),
                    Month = as.numeric(format(Date, format = "%m")),
                    Day = as.numeric(format(Date, format = "%d")))
joined <- merge(joined, pressuredata, sort = FALSE, all.x = TRUE)

# avoid ambiguity and change colname for met wind; add Date column
names(windsdata)[6] <- "metWind"

toDate <- function(year, month, day) {
  ISOdate(year, month, day)
}

windsdata$Date <- with(windsdata, toDate(Year, Month, Day))
windsdata <- transform(windsdata, Date = as.POSIXct(as.character(Date), format = "%Y-%m-%d"))


## create appropriate averages for each instance (i.e. sampling date minus 14 days)
## FIXME: tinker with number of days, current selection give 13 previous days 
## hugely complicated scenario where appropriate values are grabbed from Lake-specific data
if(regina) {
  windsub <- subset(windsdata, Superstation == "regina")
  windgrab <- function(interval, regina) {
    datalist <- list()
    datalist <- windsub$metWind[windsub$Date < interval[[2]][1] & windsub$Date > interval[[2]][2]]
    }
} else {
  print("to be developed") # i.e. case where, despite grabbing pressure data from regina, grabbing 
  #   wind data from individual weather stations
}

if(regina) {
  windsub <- subset(windsdata, Superstation == "regina")
  spldf <- with(joined, split(joined, list(Lake), drop = TRUE)) # split into list by lakes
  sampledatelist <- lapply(spldf, "[", c("Lake", "Date")) # grab date and lake from each lake
  windraw <- list()
  windmeans <- list()
  winddf <- data.frame(Lake = character(0), sampleDate = character(0), meanWind = integer(0))
  for (i in 1:length(sampledatelist)) {
      sampledates <- sampledatelist[[i]] 
      interval <- list()
          for (i in 1:nrow(sampledates)) { 
          interval[[i]] <- list(sampledates[['Lake']][1], seq.POSIXt(from = sampledates[i,'Date'], 
                                                                     by = "-14 day", length.out = 2))
          } # this gives a list the length of nrow(sampledates[[i]]) with lake id and date range
      windraw <- lapply(interval, windgrab) # gives list length of nrow(sampledates[[i]]) 
      #   with all 13 wind data (with setting -14 days)
      windmeans <- unlist(lapply(windraw, mean, na.rm = TRUE))
      ## FIXME: is this how function options are entered in lapply?
      winddf <- rbind(winddf, data.frame(Lake = sampledates['Date'], sampleDate = sampledates['Lake'], 
                                         meanWind = windmeans)) 
  }
} else {
  print("not-regina to be developed") ## FIXME not-regina option null
}

joinedtest <- merge(joined, winddf, by = c("Date", "Lake"), sort = FALSE, all.x = TRUE)
# one data frame with means for all sample occasions of all lakes
## FIXME: If joinedtest is deemed fully operational, the code below needs to be changed and run to
##    produce a new gasFlux.rds, also to incorporate the change in wind transformation (below).for now,
##    joinedtest will be saved as a separate object which only contains parameters prior to function 
##    executions (since this can then be run separately under different scenarios in co2_scenarios.R)

## need to change wind from km/h to m/s
## NOTE! chatted with kerri and it seems the database is wrong! despite column name km/h the data 
##    are actually m/s --> commented out Wind conversion below.
# joined <- transform(joined, Wind = Wind * (1000/60/60))
joinedtest <- transform(joinedtest, meanWind = meanWind * (1000/60/60))

## need to change dic from mg/L to uM
joined <- transform(joined, TIC = TIC / 0.012)
joinedtest <- transform(joinedtest, TIC = TIC / 0.012)
#   this is what Kerris' spreadsheet indicates for the unit conversion

## save joinedtest for now
saveRDS(joinedtest, "data/private/joinedtest.rds")

## need to source function that will run through the calculations
source("functions/gasExchange.R") ## FIXMEs for scenario = kerri

## Run the gas exchange equations on our data to create co2 flux
joined <- transform(joined,
                    co2Flux = gasExchange(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC, 
                                          pco2atm = rep(370, times = 1242), # see line 94
                                          kpa = Pressure, wind = Wind, salt = Salinity))
#     archaic issue with naming in the database, dic = TIC is correct

## Save output for paper
saveRDS(joined, "data/private/gasFlux.rds")
