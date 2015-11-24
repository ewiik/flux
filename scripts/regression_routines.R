## regression model development for routines CO2 flux

## get necessary explanatory data sets
co2expl <- readRDS("data/private/co2explained.rds") # YEAR, Month, LAKE, each sampling date
precip <- readRDS("data/precip.rds") # Year, Month, for snow and rain and both
nao <- readRDS("data/naoseasonal.rds") # trimonth subsets for a whole year; all colnames lowercase
#   nao not in model right now
soistand <- readRDS("data/soi_stand.rds") # monthly; all colnames uppercase
pdo <- readRDS("data/pdo.rds") # monthly, and annual mean: all colnames lowercase

## get CO2 flux estimates
fluxes <- readRDS("")

## choose params that we want in model from co2expl offerings
## at the mo, relative humidity is in as a proxy for evaporation effects.. Considered this better than using
##    just precipitation since that may not be directly linked due to advective-dominated precip
regvars <- subset(co2expl, select = c("YEAR", "Month", "Date","DOY", "LAKE", "Chl_a_ug_L", "lakeNPP",
                                      "lakeR", "TDN_ug_L", "pH_surface", "DOC_mg_L", "Oxygen_ppm",
                                      "AirTempMonthly", "RelHum"))

## change colnames and designations to be compatible
colnames(regvars)[which(colnames(regvars) == "YEAR")] <- "Year"
colnames(pdo)[which(colnames(pdo) == "year")] <- "Year"
pdostack <- stack(pdo, select = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", 
                                  "sep", "oct", "nov", "dec"))
