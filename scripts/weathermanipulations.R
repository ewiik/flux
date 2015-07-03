## we have 1080 dfs in total; in met, they do not have the first 9 lines
## grab met from gasexchange_pressuredata.R

library("reshape")

## get station details
stations <- read.csv("data/weatherstations.csv")

## here something to add column "station"
names(met[[1]])
# [1] "StationID"             "Date/Time"             "Year"                  "Month"                
# [5] "Day"                   "Time"                  "Data Quality"          "Temp (degC)"          
# [9] "Temp Flag"             "Dew Point Temp (degC)" "Dew Point Temp Flag"   "Rel Hum (%)"          
# [13] "Rel Hum Flag"          "Wind Dir (10s deg)"    "Wind Dir Flag"         "Wind Spd (km/h)"      
# [17] "Wind Spd Flag"         "Visibility (km)"       "Visibility Flag"       "Stn Press (kPa)"      
# [21] "Stn Press Flag"        "Hmdx"                  "Hmdx Flag"             "Wind Chill"           
# [25] "Wind Chill Flag"       "Weather"  

# met2 <- met[c(1:2, 1000:1002)] # testset, has regina and indhead
metmelt <- melt.list(met, id.vars = c("StationID", "Date/Time"))

metmelt$superstation[metmelt$StationID == 3002 | metmelt$StationID == 51441] <- "regina"
metmelt$superstation[metmelt$StationID == 3062 | metmelt$StationID == 44203 | metmelt$StationID == 48977] <- "yorkton"
metmelt$superstation[metmelt$StationID == 3318] <- "outlook"
metmelt$superstation[metmelt$StationID == 2926 | metmelt$StationID == 2925] <- "indhead"
unique(metmelt$superstation) # yes all appear, no NAs

metsub <- subset(metmelt, subset = grepl("kPa", variable), select = c(variable,value)) # works
metsubregina <- subset(metmelt, subset = c(superstation == "regina", grepl("kPa", variable)), select = c(variable,value)) # doesn't work
metsubregina <- subset(metmelt, subset = c(superstation == "regina", variable = grepl("kPa", variable)), 
                       select = c(variable,value)) # crashes computer

  