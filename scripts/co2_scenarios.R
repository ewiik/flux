## received complete time series of data from kerri (.../fromkerri/): pco2foremmasept2015 
##    and pco2foremmatidy - the latter has been used for sheet1, and sheet 2 grabbed from former,
##    saved as "data/private/kerri_co2flux_results_complete.csv" and 
##    "data/private/kerri_co2flux_results_complete_2.csv"
## these files have some assumptions that differ from the original calcs in the Nature pub;
##    this script intends to reproduce exactly what was done before (rather than check
##    reproducibility of code as per co2data_comparisons.R) and also allow for comparison of
##    conclusions of regressions run with different options. Therefore necessary columns have been added
##    to params (wind and altitude) or embedded in gasExchangeFlex.R (salt = 0, missing DIC, pco2atm)
## altitude values from /fromkerri 2007 dic based on cond xls are identical to Elevation in params
## FIXME: still need to sort out removing outliers consistently across scripts (incl. conductivity, salinity)

## NOTE: *DIC values*: Kerri did DIC (mg/L) = 26.57 + 0.018 * Conductivity (see email 1st June 2015)
##    for missing DIC data points! 
##       *CO2uM, pCO2*: kerri used the spreadsheet to calculate from DIC, but when DIC missing 
##    (even with the conductivity relationship), she used the final regression of pH to CO2 and pH to pCO2 
##    in order to fill in the missing numbers - so this should only have been evoked when DIC & cond missing:
##    log CO2 (uM) = 10.94 -1.124*pH  (r2 = 0.94)
##    log pCO2 (uatm) = 12.076-1.101*pH (r2 = 0.95)
## ergo since we using only existing data, the CO2 stuff is ok. 

## get necessary data 
ml <- read.csv("data/maunaloa.csv") 

if (file.exists("data/private/params.rds")) {
  params <- readRDS("data/private/params.rds")  
} else {
  stop("run pressuremanipulations.R first as 'params.rds' doesn't exist")
  # FIXME this should really source("pressuremanipulations.R")
}

## Create necessary values for starting parameters
# insert pco2atm from Mauna Loa
mlsub <- subset(ml, select = c('year', 'month', 'value'))
names(mlsub) <- c("Year", "Month", "pco2atm")
params <- merge(params, mlsub, by = c("Year", "Month"))

# insert kerri's wind values (finlayetal2009 p. 2556); Pasqua not done here but I calculated average 
#   as 3.7. NOTE: if I calculate means of available data I have for all lakes, the values are different.
winds <- data.frame(Lake = c("B", "C", "D", "K", "L", "P", "WW"),
                    Wind_ms = c(4.1, 4.2, 3.4, 3.3, 4.3, 3.7, 2.8)) 
params <- merge(params, winds, by = "Lake") # adds Kerri's wind values
replace <- which(names(params) %in% c("Wind", "Wind_ms"))
names(params)[replace] <- c("measuredWind", "kerriWind") # first one is the wind of the *sampling date*      

## remove true outliers where pH > 12
params <- subset(params, pH < 12)

## source function that will run through the calculations
source("functions/gasExchangeFlex.R") 

## create our two scenarios; parameters in function = temp, cond, ph, wind, kerri = FALSE, salt = NULL, 
##    dic = NULL, alt = NULL, kpa = NULL, pco2atm = NULL, trace = FALSE
## FIXME: test more scenarios, double check code, tweak.
## CHECKED: ran co2data_comparisons.R, took joined object, plotted co2Flux from joined and params 
##    (scenario = new) and the data are nearly reproducing!! (with the archaic flux rds, i.e. different
##    wind data)
scenario <- "new"
co2Flux <- switch(scenario,
                  new      = with(params, gasExchangeFlex(Temperature, Conductivity, pH, meanWind, kerri= FALSE, 
                                             salt = Salinity, dic = TIC, kpa = Pressure, 
                                             pco2atm = pco2atm)),
                  original = with(params, gasExchangeFlex(Temperature, Conductivity, pH, kerriWind, kerri= TRUE, 
                                             dic = TIC, kpa = Pressure, 
                                             alt = Elevation)))
params$co2Flux <- co2Flux
