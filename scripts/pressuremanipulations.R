## outputting co2 exchange with all data combined as one table
## get necessary data from database queries

if (file.exists("data/pressuredata.rds")) {
  pressuredata <- readRDS("data/pressuredata.rds")  
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

## merge data on weather superstation for each lake
joined <- merge(joined, superstations, sort = FALSE)

## merge joined with pressure data
## add Year & Month fields to aid the match
##
## FIXME: Still need to decide what to do with missing pressure data.
##        Do we use pressure data from other near site? If so, why not use
##        that for all the data not just missing stuff?
joined <- transform(joined, Year = as.numeric(format(Date, format = "%Y")),
                    Month = as.numeric(format(Date, format = "%m")))
joined <- merge(joined, pressuredata, sort = FALSE, all.x = TRUE)

## need to change wind from km/h to m/s
joined <- transform(joined, Wind = Wind * (1000/60/60))

## need to change dic from mg/L to uM
joined <- transform(joined, TIC = TIC / 0.012)
#   this is what Kerris' spreadsheet indicates for the unit conversion

## need to source function that will run through the calculations
source("functions/gasExchange.R") #FIXMEs in here too!!!

## Run the gas exchange equations on our data to create co2 flux
joined <- transform(joined,
                    co2Flux = gasExchange(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC, 
                                          pco2atm = rep(370, times = 1242), # see line 94
                                          kpa = Pressure, wind = Wind, salt = Salinity))
#     archaic issue with naming in the database, dic = TIC is correct

##      FIXME: Still need to decide what to do with pco2atm - Do we use blanket 370 like Kerri did? Or
## ftp://aftp.cmdl.noaa.gov/data/trace_gases/co2/in-situ/surface/mlo/co2_mlo_surface-insitu_1_ccgg_MonthlyData.txt
##      i.e. data/maunaloa.csv

## For this work we only want a subset
take <- c("B", "C", "K", "L", "P", "WW")
gasFlux <- subset(joined, subset = Lake %in% take)

## Save output for paper
saveRDS(gasFlux, "data/private/gasFlux.rds")
