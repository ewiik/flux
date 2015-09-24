## received complete time series of data from kerri (.../fromkerri/): pco2foremmasept2015 
##    and pco2foremmatidy - the latter has been used for sheet1, and sheet 2 grabbed from former
## these files have some assumptions that differ from the original calcs in the Nature pub;
##    this script intends to reproduce exactly what was done before (rather than check
##    reproducibility of code as per co2data_comparisons.R) and also allow for comparison of
##    conclusions of regressions run with different options.
## this script will require gasExchangeFlex.R
## FIXME: still need to sort out removing outliers consistently across scripts

## get necessary data from previous processing
## FIXME: will this be necesssary?
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

## NOTE: *DIC values*: Kerri may have done DIC (mg/L) = 26.57 + 0.018 * Conductivity (see email 1st June 2015)
##    for missing, or all, DIC data points! she is resolving this (see email 23.09.2015)
##       *CO2uM, pCO2*: kerri used the spreadsheet to calculate from DIC, but when DIC missing 
##    (even with the conductivity relationship), she used the final regression of pH to CO2 and pH to pCO2 
##    in order to fill in the missing numbers - so this should only have been evoked when DIC & cond missing:
##    log CO2 (uM) = 10.94 -1.124*pH  (r2 = 0.94)
##    log pCO2 (uatm) = 12.076-1.101*pH (r2 = 0.95)
## ergo since we using only existing data, the CO2 stuff is ok. However: 
## FIXME: Kerri to report back on DIC values

## Grab starting parameters and create a few more necessary columns
params <- readRDS("data/private/joinedtest.rds") 
ml <- read.csv("data/maunaloa.csv") 
mlsub <- subset(ml, select = c('year', 'month', 'value'))
names(mlsub) <- c("Year", "Month", "pco2atm")
params <- merge(params, mlsub, by = c("Year", "Month"))

## need to source function that will run through the calculations
source("functions/gasExchangeFlex.R") 

## test whether the function is working; params = temp, cond, ph, wind, kerri = FALSE, salt = NULL, dic = NULL, 
##    alt = NULL, kpa = NULL, pco2atm = NULL
## FIXME: test more scenarios, double check code, tweak.
## scenario where all should work
flexing <- transform(params, co2Flux = gasExchangeFlex(Temperature, Conductivity, pH, meanWind, kerri= FALSE, 
                                                       salt = Salinity, dic = TIC, kpa = Pressure, 
                                                       pco2atm = pco2atm))
## kerri = TRUE, failing to provide alt
flexing <- transform(params, co2Flux = gasExchangeFlex(Temperature, Conductivity, pH, Wind, kerri= TRUE, 
                                                       salt = Salinity, dic = TIC, kpa = Pressure, 
                                                       pco2atm = pco2atm))
## kerri = TRUE, provide alt
flexing <- transform(params, co2Flux = gasExchangeFlex(Temperature, Conductivity, pH, Wind, kerri= TRUE, 
                                                       salt = Salinity, dic = TIC, kpa = Pressure, 
                                                       alt = Elevation))

## kerri = TRUE, no salt provided
flexing <- transform(params, co2Flux = gasExchangeFlex(Temperature, Conductivity, pH, Wind, kerri= TRUE, 
                                                       dic = TIC, kpa = Pressure, 
                                                       alt = Elevation))

