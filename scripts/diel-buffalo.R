## load in Helen's buoy data -- from file sent 28.01.2016
## two sonde depths, n.1 is the shallower one based on PAR..
## Not sure what cdom and mv are measures of.. and why so many temps?
## FIXME: really need to learn what's what with the column name matching
## FIXME: removed outliers but also 0s that are not real 0s.. need to remove

## load necessary packages
library("ggplot2")

## read in data
bdat <- read.csv("../data/private/BPBuoyData2015raw.csv", skip = 4, 
                 col.names = c("datetime", "batvolt", "winddir", "windsp", "airtemp", 
                               "relhum", "pressure", "dailyrain", "par1", "par2",
                               "par3", "cdom", "mvolts", "co2-1", "co2-2", "temp1",
                               "cond1", "ph1", "ph1mV", "turb", "chl", "chlrfu", "bga1cell",
                               "bga1rfu", "ODOrel1", "ODOabs1", "temp2", "cond2", "ph2", 
                               "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3", 
                               "temp4", "temp5", "temp6", "temp7"))

## make time understandable
bdat <- transform(bdat, datetime = as.POSIXct(as.character(datetime), 
                                              format = "%m/%d/%Y %I:%M:%S %p"))
bdat <- transform(bdat, Hour = as.numeric(format(datetime, format = "%H")))
bdat <- transform(bdat, Month = as.numeric(format(datetime, format = "%m")))
bdat <- transform(bdat, Day = ifelse(Hour < 21 & Hour > 7, TRUE, FALSE))

## check the CO2 columns
with(bdat, plot(co2.1 ~ datetime))

## not sure why all june and july values are completely different for both co2 profiles..!?     
with(bdat, plot(density(co2.1)))
with(bdat, plot(density(co2.2)))

## remove the crazy negative and positive values and bound the ppm
##    but first make note of which rows these are
outliers <- which(abs(bdat$co2.1) > 1200 | abs(bdat$co2.2) > 1200)
bdat<- bdat[-outliers,]
## still bimodal but don't wanna remove too much in case it's real

## plot some diagnostics and look if pH and oxygen determining CO2
##    assuming I know which odo and ph measure relates to which co2 probe...
with(bdat, plot(co2.1 ~ Hour, col = ifelse(Day, "red", "black")))
with(bdat, plot(co2.2 ~ Hour, col = ifelse(Day, "red", "black")))

with(bdat, plot(co2.1 ~ ph1, col = ifelse(Month < 9, "red", "black")))
with(bdat, plot(co2.1 ~ ODOrel1, col = ifelse(Day, "red", "black")))
## something wrong with ODO measurements!! loads 0s. Let's remove 0s
## Existing data mostly contiguous in time
with(bdat[-which(bdat$ODOrel1 == 0), ], plot(co2.1 ~ ODOrel1, col = ifelse(Day, "red", "black")))


with(bdat, plot(co2.2 ~ ph2, col = ifelse(Month < 9, "red", "black")))
with(bdat, plot(co2.2 ~ ODOrel2, col = ifelse(Month < 9, "red", "black")))
## here, something fishy with the month of September.. CO2 mostly 0

with(bdat[-which(bdat$ODOrel1 == 0), ], 
     plot(ODOrel1 ~ ph1, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
with(bdat, plot(ODOrel2 ~ ph2, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
## night and day different slopes in sonde1, doesn't happen in sonde2 even when
##    following the ODOrel1 subset

with(bdat[-which(bdat$ODOrel1 == 0), ], 
     plot(ODOrel1 ~ chl, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
with(bdat, plot(ph1 ~ chl, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
