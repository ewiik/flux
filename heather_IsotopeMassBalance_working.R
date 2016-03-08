# Source FUnctions 
source("C:/Users/limno/Dropbox/daily transfers between mac and PC/Chapter 1/R Code/isotopeFun_Jan6_H.R")

#import Residence time data
residencetime <- read.csv("C:/Users/limno/Dropbox/daily transfers between mac and PC/Chapter 1/Data/residencetime.csv")
#View(residencetime)

#import water isotope data
qisodata <- read.csv("C:/Users/limno/Dropbox/daily transfers between mac and PC/Chapter 1/Data/qisodata.csv")
#View(qisodata)

#import meterological data- humidity and temperature
metdata <- read.csv("C:/Users/limno/Dropbox/daily transfers between mac and PC/Chapter 1/Data/metdata2.csv")
#View(metdata)

#import precipitation data
metprecip <- read.csv("C:/Users/limno/Dropbox/daily transfers between mac and PC/Chapter 1/Data/precipquappelle95_2014csv.csv")
#View(precipquappelle95_2014csv)
metprecip <- transform(metprecip, Date.Time = as.Date(as.character(Date.Time), format = "%m/%d/%Y"))
names(metprecip)[6:8] <- paste("Total", c("Rain","Snow","Precip"), sep = ".")
metprecip <- transform(metprecip, StationID = gsub(" $", "", as.character(StationID)))

#import lummination data
lummination <- read.csv("C:/Users/limno/Dropbox/daily transfers between mac and PC/Chapter 1/Data/lummination.csv")
#View(lummination)

#import climate stations
stationID <- read.csv("C:/Users/limno/Dropbox/daily transfers between mac and PC/Chapter 1/Data/stationID.csv")
#View(stationID)

#Change date for delL to %m-%d-%Y using as.date and POISIX ct NOT WORKING... No idea what this code is suppose to do
qisodata <- transform(qisodata,
                      sampleDate = as.Date(as.character(date),
                                           format = "%m/%d/%Y"))

## trying to convert dates but I dont know what I am doing wrong 
metdata <- 
  transform(metdata,
            Date.Time = as.POSIXct(as.character(Date.Time),
                                   format = "%Y-%m-%d %H:%M"))

## merge the station name with the met data by sation id
metdata <- merge(metdata, stationID, sort = FALSE)

## Drop all met data prior to 1967; we do this because there was DST in Sask
## up to 1966
take <- metdata$Year < 1967
metdata <- metdata[!take, ] # not take == drop

## Need a sampleDay field to aggregate Temp & RelHum by day in the loop
metdata <- transform(metdata,
                     sampleDay = as.Date(Date.Time))

## Merge residence time data with isotope data, and...
qisodata <- merge(qisodata, residencetime, sort = FALSE)
## ... create endpoint dates all at once
qisodata <- transform(qisodata,
                      endpoint = sampleDate - resdays)

## Monthly lummination
lummin <- aggregate(Total.hours.of.illumination ~ Month, data = lummination, FUN = mean)
names(lummin)[2] <- "illumin"

## yearly heat index
metdata <- transform(metdata, Year = as.numeric(Year),
                     Month = as.numeric(Month), Day = as.numeric(Day))
monthlyTemp <- aggregate(Temperature ~ StationName + Year + Month,
                         data = metdata, FUN = mean)
monthlyTemp <- transform(monthlyTemp, Talpha = Temperature,
                         MonthC = month.abb[Month], stringsAsFactors = FALSE)
monthlyTemp$Talpha[monthlyTemp$Talpha < 0] <- 0
## Merge with lummination data
monthlyTemp <- merge(monthlyTemp, lummin, by.x = "MonthC", by.y = "Month", sort = FALSE)

## Number of days per month
daysInMonth <- c(31,28.5,31,30,31,30,31,31,30,31,30,31)
names(daysInMonth) <- month.abb
monthlyTemp <- transform(monthlyTemp,
                         daysInMonth = daysInMonth[as.character(MonthC)])

## Iheat for each month
heatIndex <- aggregate(Talpha ~ StationName + Year, data = monthlyTemp,
                       FUN = Iheat, strict = FALSE)
names(heatIndex)[3] <- "Iheat"

## ...now merge monthlyTemp and heatIndex
monthlyTemp <- merge(heatIndex, monthlyTemp, sort = FALSE)
## resort
monthlyTemp <- monthlyTemp[with(monthlyTemp, order(StationName, Year, Month)), ]

## Add PET variable to monthlyTemp
monthlyTemp <- transform(monthlyTemp, PET = PET(illumin, Talpha, daysInMonth, Iheat))
## This works, except for years with incomplete months - esp winter months only
## then Talpha will be 0 for those months -> iheat being 0 for that year,
## then we divide by zero hence the NaN
## This is OK though as those yaers are well before the met data we need for isotope
## samples.
## **THIS ONLY AFFECTS INDIAN HEAD in 1995 Nov & Dec**

## Prepare isotiope data object with 2 columns to receive the flux weighted values
##Making the two new columns for winter and summer weather conditions
qisodata <- cbind(qisodata,
                  fluxTemp = rep(NA, nrow(qisodata)),
                  fluxRelHum = rep(NA, nrow(qisodata)),
                  winterPrecip = rep(NA, nrow(qisodata)), 
                  summerTemp = rep(NA, nrow(qisodata)))

## progressbar
pb <- txtProgressBar(min = 0, max = nrow(qisodata), style = 3)

#Loop
for (i in seq_len(nrow(qisodata))) {
  setTxtProgressBar(pb, i)
  ## lake code for current isotope sample
  curlake <- as.character(qisodata$lake[i])
  
  ## met station name for current sample
  curstation <- as.character(qisodata$StationName[i])
  
  ## sampleDate & sample time (made up) for current sample
  sampleDate <- qisodata$sampleDate[i]
  sampleTime <- as.POSIXct(paste(as.character(sampleDate), "23:59:59"))
  
  ## endpoint date & endpoint time (made up)for current sample 
  endpoint <- qisodata[i, "endpoint"]
  endtime <- as.POSIXct(paste(as.character(endpoint), "00:00:00"))
  
  ## Determine year of sample date
  sampleYear <- format(qisodata$sampleDate[i], "%Y")
  
  ## Get the year previous to the sample
  PreYear <- format(qisodata$sampleDate[i]-365, "%Y")
  
  ## Make a variable that will hold all the precip data from the desiered time period
  ## Dec Previous Year to March sample year
  
  winterStart <- as.Date(paste(PreYear, "12", "01", sep = "-"))
  winterEnd <- as.Date(paste(sampleYear, "03", "31", sep = "-"))
  summerStart <- as.Date(paste(sampleYear, "05", "01", sep = "-"))
  summerEnd <- as.Date(paste(sampleYear, "08", "31", sep = "-"))
  
  winterTake <- with(metprecip, (Date.Time >= winterStart & Date.Time <= winterEnd) &
                       StationID == curstation)
  winterPrecip <- sum(metprecip[winterTake, "Total.Precip"], na.rm = TRUE)
  
  ##Make a variable that will hold all the met data from the desiered time period
  ## May 1st to Aug 31 of sample year
  
  summerTake <- with(metdata, (sampleDay >= summerStart & sampleDay <= summerEnd) &
                       StationName == curstation)
  summerTemp <- mean(metdata[summerTake, "Temperature"], na.rm = TRUE)
  
  #Find summer precip
  #summerPrecip <- with(metprecip, (Date.Time >=summerStart  & Date.Time <= summerEnd) &
  # StationID == curstation)
  
  #summerTemp <- sum(metdata[summerPrecip, "Total.Precip"], na.rm = TRUE)
  
  ## Stick the winter and summer met data into the qisodata data frame
  qisodata[i, "summerTemp"] <-  summerTemp
  qisodata[i, "winterPrecip"] <- winterPrecip
  
  ## Met data for the current sample for its antecedent period
  take <- (metdata$Date.Time <= sampleTime & 
             metdata$Date.Time >= endtime) &
    metdata$StationName == curstation
  curmet <- metdata[take, , drop = FALSE]
  
  if (isTRUE(identical(nrow(curmet), 0L)))
    next
  
  ## Compute the daily average temperature & RelHum for current met
  dailyMet <- aggregate(cbind(Temperature, RelativeHumidity) ~ sampleDay, data = curmet,
                        FUN = mean)
  dailyMet <- transform(dailyMet, Month = format(sampleDay, format = "%b"),
                        Year = format(sampleDay, format = "%Y"), stringsAsFactors = FALSE)
  
  ## Generate weighted Temperature & RelativeHumidity for each day,
  ## weights given by PET
  dailyMetWtd <- merge(dailyMet, monthlyTemp[monthlyTemp$StationName == curstation, 
                                             c("PET","Year","MonthC","daysInMonth")],
                       by.x = c("Year","Month"), by.y = c("Year", "MonthC"))
  
  ## Add wtd Temperature and Rel humidity to dailyMetWtd:
  ## Spreading PET (mm per month) equally over the number of days in each month so we
  ## can apply this daily PET to each daily met data value
  dailyMetWtd <- transform(dailyMetWtd,
                           wtdTemp = (PET / daysInMonth) * Temperature,
                           wtdRelHum = (PET / daysInMonth) * RelativeHumidity)
  
  ## Flux weighted values for this sample
  fluxTemp <- with(dailyMetWtd, sum(wtdTemp) / sum(PET / daysInMonth))
  fluxRelHum <- with(dailyMetWtd, sum(wtdRelHum) / sum(PET / daysInMonth))
  
  ## Stick those into the isotope data frame
  qisodata[i, "fluxTemp"] <- fluxTemp
  qisodata[i, "fluxRelHum"] <- fluxRelHum
  
  ## End the loop
}
close(pb)

## Convert flux Temp to Kelvin & fluxRelHumid to decimal notation
qisodata <- transform(qisodata,
                      fluxTemp  = fluxTemp + 273.15,
                      fluxRelHum = fluxRelHum / 100)

## Weighted Average Precipitation
delP18O <- -14.6519 
delP2H <- -114.456
qisodata <- transform(qisodata,
                      delP18O = rep(delP18O, nrow(qisodata)),
                      delP2H = rep(delP2H, nrow(qisodata)))

## Use Mass Balance functions to find del E and del * for both 2H and 18O for each point
mb18O <- with(qisodata, MassBalance18O(fluxRelHum, fluxTemp, delP18O, corrected18O))
qisodata <- cbind(qisodata, as.data.frame(mb18O))
mb2H <- with(qisodata, MassBalance2H(fluxRelHum, fluxTemp, delP2H, correctedD))
qisodata <- cbind(qisodata, as.data.frame(mb2H))

## Function to fit linear regression lines to triplets of isotope data
isoLreg <- function(iso) {
  y <- c(iso$delE2H, iso$correctedD, iso$delStar2H)
  X <- cbind(rep(1,3), c(iso$delE18O, iso$corrected18O, iso$delStar18O))
  fit <- lm.fit(x = X, y = y)
  betas <- unname(coef(fit))
  list(intercept = betas[1], slope = betas[2])
}

## Use our regression function on each row - need a wrapper which is the anonymous function
coefs <- lapply(seq_len(nrow(qisodata)), function(i, data) { isoLreg(data[i, , drop = FALSE]) },
                data = qisodata)
coefs <- matrix(unlist(coefs), ncol = 2, byrow = TRUE)
colnames(coefs) <- c("intercept", "slope")

## the intercept and slope data to the isotope data
qisodata <- cbind(qisodata, coefs)

isoInputW <- function(iso) {
  isointercept(iso$slope, iso$intercept)
}

inputs <- sapply(seq_len(nrow(qisodata)), function(i, data) { isoInputW(data[i, , drop = FALSE]) },
                 data = qisodata)
inputs <- t(inputs)
colnames(inputs) <- c("delI18O", "delI2H")
qisodata <- cbind(qisodata, inputs)

#calculate E:I 
qisodata <- transform(qisodata,
                      EI18O = EI18O(fluxRelHum, ek18O, eps18O, corrected18O, delI18O, delStar18O),
                      EI2H = EI2H(fluxRelHum, ek2H, eps2H, correctedD, delI2H, delStar2H))

## Save this as a serialized object
saveRDS(qisodata, file = "quappelle-processed-isodata.rds")
