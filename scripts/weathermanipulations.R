## we have 1080 dfs in total; in met, they do not have the first 9 lines
## grab met from gasexchange_pressuredata.R

## here something to add column "station"
names(met[[1]])
# [1] "StationID"             "Date/Time"             "Year"                  "Month"                
# [5] "Day"                   "Time"                  "Data Quality"          "Temp (degC)"          
# [9] "Temp Flag"             "Dew Point Temp (degC)" "Dew Point Temp Flag"   "Rel Hum (%)"          
# [13] "Rel Hum Flag"          "Wind Dir (10s deg)"    "Wind Dir Flag"         "Wind Spd (km/h)"      
# [17] "Wind Spd Flag"         "Visibility (km)"       "Visibility Flag"       "Stn Press (kPa)"      
# [21] "Stn Press Flag"        "Hmdx"                  "Hmdx Flag"             "Wind Chill"           
# [25] "Wind Chill Flag"       "Weather"  

met2 <- met
sapply(met2[[i]], met2[[i]]$superstation = ))