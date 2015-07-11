## today we will begin work on outputting co2 exchange
## get necessary data from database queries

if (file.exists("pressuredata.rds")) {
  pressuredata <- readRDS("pressuredata.rds")  
} else {
  print("run weathermanipulations.R")
  }

if (file.exists("data/fluxquery1.csv")) { # DIC, pH, wind
  surfd <- read.csv("data/fluxquery1.csv")
} else {
  print("get fluxquery1.csv from dropbox")
}

if (file.exists("data/fluxquery2.csv")) { # elevation and lake 
  # abbreviation data
  elevd <- read.csv("data/fluxquery2.csv")
} else {
  print("get fluxquery2.csv from dropbox")
}

if (file.exists("data/fluxquery3.csv")) { # T, cond, salinity
  tcsd <- read.csv("data/fluxquery3.csv")
} else {
  print("get fluxquery3.csv from dropbox")
}

## are we really removing all NA data? THERE'S LOTS!
miss <- is.na(surfd)
miss <- rowSums(miss) != 0
surfset <- droplevels(surfd[!miss, ]) 

miss <- is.na(tcsd)
miss <- rowSums(miss) != 0
tcsset <- droplevels(tcsd[!miss, ]) 

## create date objects for later?
tcsset$datect <- as.POSIXct(as.character(tcset$Date), format = "%d-%m-%Y")

surfset$datect <- as.POSIXct(as.character(surfset$Date), format = "%d-%m-%Y")

## need to match all individual observations with correct lake and date
## still probs very wrong but will figure it out
matchtcs <- subset(tcsset, LAKE == "B" & Date %in% surfset$Date)
rownames(matchtcs) <- matchtcs$Date
matchsurf <- subset(surfset, LAKE == "B" & Date %in% tcsset$Date)
rownames(matchsurf) <- matchsurf$Date
matched <- cbind(matchtcs, matchsurf)

## need to change wind from km/h to m/s

## need to generate function that will run through the stuff in gasexchange_equations.R
source("functions/gasExchange.R")

