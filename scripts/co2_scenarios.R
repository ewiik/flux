## received complete time series of data from kerri (.../fromkerri/): pco2foremmasept2015 
##    and pco2foremmatidy - the latter has been used for sheet1, and sheet 2 grabbed from former
## these files have some assumptions that differ from the original calcs in the Nature pub;
##    this script intends to reproduce exactly what was done before (rather than check
##    reproducibility of code as per co2data_comparisons.R) and also allow for comparison of
##    conclusions of regressions run with different options.
## FIXME: meet with Kerri to discuss exact replication (see co2data_comparisons.R and gasExchange.R)

## get necessary data from previous processing
if (file.exists("data/private/gasFlux.rds")) {
  co2emma <- readRDS("data/private/gasFlux.rds")  
} else {
  print("run pressuremanipulations.R")
}

## recreate old, i.e. grab complete param data from Excel sheets and change params that differ
co2kerri <- read.csv("data/private/kerri_co2flux_results_complete.csv") # most params
names(co2kerri) <- c("Lake", "Date", "pH", "Alkalinity", "Temperature", "Conductivity", "pCO2atm", "Pressure", 
                     "Altitude", "Ionstrength", "pk1", "pk2", "ao", "a1", "a2", "pkh", "DIC", "CO2uM", "CO2uatm", 
                     "CO2eq", "CO2log", "CO2sat", "alpha", "Flux", "FluxEnh")
co2kerri <- transform(co2kerri, Date = as.POSIXct(as.character(Date), format = "%d-%m-%Y"))
co2kerri <- transform(co2kerri, Year = as.numeric(format(Date, format = "%Y")),
                      Month = as.numeric(format(Date, format = "%m")))

co2kerri2 <- read.csv("data/private/kerri_co2flux_results_complete_2.csv") # salinity, wind, temp
co2kerri2 <- transform(co2kerri2, Date = as.POSIXct(as.character(Date), format = "%d-%m-%Y"))
co2kerri <- merge(co2kerri, co2kerri2, by = c("Lake", "Date"), all.x = TRUE)

# drop calculated variables
co2kerri <- subset(co2kerri, select = c('Lake', 'Date', 'Year', 'Month', 'pH', 'Temperature'
                                        , 'Conductivity', 'pCO2atm', 'Altitude', 'DIC', 'Wind_ms'
                                        , 'Salinity_ppt')) # recall Pressure not originally used,
                                        # and Alkalinity calculated from DIC  

# remove true outliers where pH > 12
co2kerri <- subset(co2kerri, pH < 12)

# remove lakes we don't care about
co2kerri <- subset(co2kerri, Lake %in% c("B", "C", "WW", "D", "K", "L", "P"))

# insert previously used wind values
winds <- data.frame(Lake = c("B", "C", "D", "K", "L", "P", "WW"), Wind_ms = c(4.1, 4.2, 3.4, 3.3, 4.3, 
                                                                              3.7, 2.8)) 
# see finlayetal2009 p. 2556; Pasqua not done here but me calculated average as 3.7. note also
#   that if I calculate means of the available data points I have for all lakes, the values are different.

# insert previously used altitude values (grabbed from /fromkerri 2007 dic based on cond xls)
altitudes <- data.frame(Lake = c("B", "C", "D", "K", "L", "P", "WW"), Altitude = c(509.3, 451.7, 552,
                                                                                   478.2, 490.1, 479, 570.5))
# merge and replace
co2kerriw <- merge(altitudes, co2kerri, by = "Lake")
co2kerri$Altitude <- co2kerriw$Altitude.x

co2kerrix <- merge(winds, co2kerri, by = "Lake")
co2kerri$Wind_ms <- co2kerrix$Wind_ms.x

# CHECK: CO2uM calculation, pCO2 calculation (one of the fixmes in gasExchange.R)
# and then by gods this data set is complete for scenario = kerri

## Grab my flux calculations and met data and start scenarios which are:
##    1: Grab pCO2 from Mauna Loa, grab wind from meteo stations and create averages for 
##    sample time - 2 weeks, grab pressure from regina which was longest-running, do same average
## Though average of time x could be a condition for the function call
## These options will appear in gasExchangeFlex.R

original <- readRDS("data/private/gasFlux.rds") # (pco2atm = 370) (see pressuremanipulations.R)
joined <- merge(co2kerrisub, original, by = c("Lake", "Date"), all.x = TRUE)

## need to source function that will run through the calculations
source("functions/gasExchange.R") 

## start messing around with parameters
## ====================================
##    Part 1: create complete similarity
##    created new gasExchangealt.R to use alt rather than kpa and suppress salt calculation
source("functions/gasExchangealt.R")
identical <- transform(joined,
                       co2Flux = gasExchangealt(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC, 
                                                pco2atm = rep(370, times = nrow(joined)), 
                                                alt = rep(509.3, times = nrow(joined)),
                                                wind = rep(4.06, times = nrow(joined)), 
                                                salt = rep(0, times = nrow(joined))))
with(identical, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(identical, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)"))
# assuming R rounding errors here

##    Part 2: introduce salt (all else same)
##    created new gasexchangesalt.R to phase back the salt calculation and use salt
source("functions/gasExchangesalt.R")
salty <- transform(joined,
                   co2Flux = gasExchangesalt(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC, 
                                             pco2atm = rep(370, times = nrow(joined)), 
                                             alt = rep(509.3, times = nrow(joined)),
                                             wind = rep(4.06, times = nrow(joined)), 
                                             salt = Salinity))
with(salty, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(salty, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)", col = Lake))
with(salty, which(FluxEnh - co2Flux == max(FluxEnh - co2Flux, na.rm = TRUE)))
# one "outlier" here is a normal day in 2001, the other one of the high pH days;
#     salt is used in equations r1 and r2
take <- with(salty, which(abs(FluxEnh - co2Flux) >= 10))
salty[take,]

##    Part 3: introduce wind (all else same)
windy <- transform(joined,
                   co2Flux = gasExchangealt(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC, 
                                            pco2atm = rep(370, times = nrow(joined)), 
                                            alt = rep(509.3, times = nrow(joined)),
                                            wind = Wind, 
                                            salt = rep(0, times = nrow(joined))))
with(windy, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(windy, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)", col = Lake))
with(windy, which(FluxEnh - co2Flux == max(FluxEnh - co2Flux, na.rm = TRUE)))
# again max difference at the pH anomaly. funnily #303 seems to produce most differences, #850 not so

##    Part 4: replace alt with Pressure (all else same)
pressed <- transform(joined,
                     co2Flux = gasExchange(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC, 
                                           pco2atm = rep(370, times = nrow(joined)), 
                                           kpa = Pressure,
                                           wind = rep(4.06, times = nrow(joined)), 
                                           salt = rep(0, times = nrow(joined))))
with(pressed, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(pressed, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)"))
# using pressure rather than altitude incurs the least changes; less than 10 units of flux.

## Save some output for later
## saveRDS(, "data/private/.rds")

