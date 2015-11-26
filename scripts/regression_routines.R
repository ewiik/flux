## regression model development for routines CO2 flux
## FIXME: do we want precip even when we have humidity? ALSO: still need to decide on what to replace 
##    ice-off with... Do we want March precip/snow AND temperature and do another interaction term?

## get necessary explanatory data sets
co2expl <- readRDS("data/private/co2explained.rds") # YEAR, Month, LAKE, each sampling date
precip <- readRDS("data/precip.rds") # Year, Month, for snow and rain and both
nao <- readRDS("data/naoseasonal.rds") # trimonth subsets for a whole year; all colnames lowercase
#   nao not in model right now
soistand <- readRDS("data/soi_stand.rds") # monthly; all colnames uppercase
pdo <- readRDS("data/pdo.rds") # monthly, and annual mean: all colnames lowercase
temp <- readRDS("data/temperaturedata.rds")# Year, Month, Superstation...

## get CO2 flux estimates (alter calculation arguments in co2_scenarios.R)
fluxes <- readRDS("data/private/params-flux.rds")

## choose params that we want in model from co2expl offerings
## at the mo, relative humidity is in as a proxy for evaporation effects.. Considered this better than using
##    just precipitation since that may not be directly linked due to advective-dominated precip
regvars <- subset(co2expl, select = c("YEAR", "Month", "Date","DOY", "LAKE", "Chl_a_ug_L", "lakeNPP",
                                      "lakeR", "TDN_ug_L", "pH_surface", "DOC_mg_L", "Oxygen_ppm",
                                      "AirTempMonthly", "RelHum"))

## do whatever needs to be done with precip data. and March things
## ============================================================================
mycolSums <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  sum <- colSums(subdf, na.rm = TRUE)
  cbind(data.frame(Year = df['Year'][1,], df['Month'][1,], t(sum)))
}

# grab total snow fall in march
snowsplit <- with(precip, split(precip, list(Year, Month), drop = TRUE))
snowtotal <- lapply(snowsplit, mycolSums, cols = c("Total Snow (cm)")) 
snowtotal <- do.call("rbind", snowtotal)
rownames(snowtotal) <- NULL
snowtotal <- subset(snowtotal, select = c(Year, Total.Snow..cm.), df..Month...1... == 3)
colnames(snowtotal)[which(colnames(snowtotal) == "Total.Snow..cm.")] <- "MarchSnowFallcm"


# grab total precip in march.. FIXME: snow-specific better?
monthsplit <- with(precip, split(precip, list(Year, Month), drop = TRUE))
monthtotal <- lapply(monthsplit, mycolSums, cols = c("Total Precip (mm)")) 
monthtotal <- do.call("rbind", monthtotal)
rownames(monthtotal) <- NULL

# grab march values for temp
marchtemp <- subset(temp, select = c(Year, Temperature), 
                                  subset = Month == 3 & Superstation == "regina")
colnames(marchtemp)[which(colnames(marchtemp) == "Temperature")] <- "MarchTemp"
## ============================================================================

## change colnames and designations to be compatible
colnames(regvars)[which(colnames(regvars) == "YEAR")] <- "Year"
colnames(regvars)[which(colnames(regvars) == "LAKE")] <- "Lake"

colnames(pdo)[which(colnames(pdo) == "year")] <- "Year"
pdostack <- stack(pdo, select = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", 
                                  "sep", "oct", "nov", "dec"))
pdostack$Year <- rep(pdo$Year, times = 12)
colnames(pdostack)[which(colnames(pdostack) == "ind")] <- "Month"
colnames(pdostack)[which(colnames(pdostack) == "values")] <- "PDO"
pdostack <- transform(pdostack, Month = rep(c(1:12), each = 116))

names(soistand) <- c("Year", 1:12)
soistack <- stack(soistand, select = -Year)
soistack$Year <- rep(soistand$Year, times = 12)
colnames(soistack)[which(colnames(soistack) == "ind")] <- "Month"
colnames(soistack)[which(colnames(soistack) == "values")] <- "SOI"

## check NA situation
fluxna <- which(is.na(fluxes$co2Flux))
nrow(fluxes[fluxna,]) # due to missing DIC mostly
nrow(fluxes) # 295 NA out of 1066

nanumbers <- rowSums(is.na(regvars)) 
dataloss <- regvars[nanumbers > 0,]
nrow(dataloss) # 278 out of 1094

## merge harmonised tables; start: 1094 rows in regvars
regvars <- merge(regvars, pdostack) # 1094
regvars <- merge(regvars, soistack) # 1094
regvars <- merge(regvars, marchtemp) # 1094
regvars <- merge(regvars, snowtotal) # 1094

vardiff <- regvars$Date[regvars$Date %in% fluxes$Date == FALSE]
fluxdiff <- fluxes$Date[fluxes$Date %in% regvars$Date == FALSE]
vardiff[order(vardiff)]
fluxdiff[order(fluxdiff)]
# seems that we don't have 1994 and 1995 in regvars --> AHH!!! production estimates only started in 1996
#   so if want those in, need to rethink....! also need to create merges that allow for this...
## FIXME: change co2explvar script merge options when doing produc + the others!!!
regvars <- merge(regvars, fluxes[,c('Lake', 'Date', 'co2Flux')])



## look at relationships etc.

