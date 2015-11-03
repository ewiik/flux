## running CO2 flux calculations for heather's data too
## FIXME: tidy up commentary, tidy up removal of carbonID --> run calcs

## download necessary files
sites <- read.csv("data/private/qCO2fluxSites.csv")
vars <- read.csv("data/private/qCO2fluxVars.csv")
alt <- read.csv("data/private/heathersitesalt.csv") # altitudes derived from a lat long website
carbon <- read.csv("data/private/carbon.csv")

## load necessary packages
library(reshape2)

## deal with units: change mS to uS 
vars$dataValue[vars$variableID == "conductivitySLS"] <- 
  vars$dataValue[vars$variableID == "conductivitySLS"] * 1000
vars$variableID[vars$variableID == "conductivitySLS"] <- "conductivity"

## remove comments section for data processing purposes
sites <- sites[,-which(colnames(sites) == "comments")]

## the tables I've joined have incompatible columns... need to individualise lakeName
## therefore renaming the second Little Manitou as Little Manitwo
sites$lakeName <- as.character(sites$lakeName) 
dups <- grep("Little Manitou", sites$lakeName)
sites[dups[2],'lakeName'] <- "Little Manitwo"

alt$lakeName <- as.character(alt$lakeName) 
dups <- grep("Little Manitou", alt$lakeName)
alt[dups[2],'lakeName'] <- "Little Manitwo"

vars$sampleID <- as.character(vars$sampleID) 
dups <- grep("Little Manitou", vars$sampleID)
dupstart <- length(dups)/2 + 1
dupend <- length(dups)
vars[dups[dupstart:dupend],'sampleID'] <- "Little Manitwo"
colnames(vars)[which(colnames(vars) == "sampleID")] <- "lakeName"

carbon$sampleID <- as.character(carbon$sampleID) 
dups <- grep("Little Manitou", carbon$sampleID)
dupstart <- length(dups)/2 + 1
dupend <- length(dups)
carbon[dups[dupstart:dupend],'sampleID'] <- "Little Manitwo"
colnames(carbon)[which(colnames(carbon) == "sampleID")] <- "lakeName"

## some data are long format, but with different lengths (nrow) so can't unstack (at least not
##    with my level of R. So, splitting by variableID to be merged....
varsplit <- with(vars, split(vars, variableID))
carbonsplit <- with(carbon, split(carbon, type))

## create function that makes each variable into a standalone column
varnames <- function(dflist, rename, renameto, remove) {
    for (i in seq(dflist)) {
      names(dflist[[i]])[which(colnames(dflist[[i]]) == rename)] <- renameto[i]
      dflist[[i]] <- dflist[[i]][-which(colnames(dflist[[i]]) %in% remove)] # == 
    }
dflist
}

## run function on varsplit, merge with sites and keep merging
names <- c(names(varsplit))
varlist <- varnames(varsplit, "dataValue", names, c("variableID", "depth"))
allvars <- sites
for (i in seq(varlist)) {
  allvars <- merge(allvars, varlist[[i]], by = "lakeName", all = TRUE)
}

## run function on carbonsplit, merge with allvars and keep merging
names <- c(names(carbonsplit))
varlist <- varnames(carbonsplit, "dataValue", names, c("type", "carbonID"))

for (i in seq(varlist)) {
  allvars <- merge(allvars, varlist[[i]], by = "lakeName", all = TRUE)
}

## merge with altitude information
allvars <- merge(allvars, alt[,c('lakeName', 'Altitude')])

## remove archaic conductivitySLS
allvars <- allvars[,-which(colnames(allvars) == "conductivitySLS")]

## need to change wind from km/h to m/s
allvars <- transform(allvars, wind = wind * (1000/60/60))

## how many NA? (first select only those vars that strictly needed for running gasflux)
dataloss <- subset(allvars, select = c(lakeName, wind, conductivity, Altitude, 
                                       pH, salinity, temperature, DIC))
nanumbers <- rowSums(is.na(dataloss))
datalost <- dataloss[rowSums(is.na(dataloss)) > 0,]
nrow(datalost) # 19; for these sites limno largely not done; no pH, mostly no temp, sal, cond

## save full dataset for running gasflux on
saveRDS(allvars, "data/private/heathergasfluxsupp.rds")

