## running CO2 flux calculations for heather's data too

## download necessary files
sites <- read.csv("../data/private/qCO2fluxSites.csv")
vars <- read.csv("../data/private/qCO2fluxVars.csv")
alt <- read.csv("../data/private/heathersitesalt.csv") # altitudes derived from a lat long website
carbon <- read.csv("../data/private/carbon.csv")
vars0 <- read.csv("../data/private/qpH@depth0.csv") # some lakes have only 0 depth pH.. need to
#   be incorporated

## load necessary packages
library(reshape2)

## deal with units: change mS to uS for SLS 
vars$dataValue[vars$variableID == "conductivitySLS"] <- 
  vars$dataValue[vars$variableID == "conductivitySLS"] * 1000
vars$variableID[vars$variableID == "conductivitySLS"] <- "conductivity"

vars0$dataValue[vars0$variableID == "conductivitySLS"] <- 
  vars0$dataValue[vars0$variableID == "conductivitySLS"] * 1000
vars0$variableID[vars0$variableID == "conductivitySLS"] <- "conductivity"

## remove comments section for data processing purposes
sites <- sites[,-which(colnames(sites) == "comments")]

## the tables I've joined have incompatible columns... need to individualise lakeName
## until database is fixed
## therefore renaming the second Little Manitou as Little Manitwo
##############################################
##############################################
##                                          ##
## FIX ME !!!!!!!!!!!!!!!!!!!!!!!           ##
##                                          ##
##############################################
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

vars0$sampleID <- as.character(vars0$sampleID) 
dups <- grep("Little Manitou", vars0$sampleID)
dupstart <- length(dups)/2 + 1
dupend <- length(dups)
vars0[dups[dupstart:dupend],'sampleID'] <- "Little Manitwo"
colnames(vars0)[which(colnames(vars0) == "sampleID")] <- "lakeName"

carbon$sampleID <- as.character(carbon$sampleID) 
dups <- grep("Little Manitou", carbon$sampleID)
dupstart <- length(dups)/2 + 1
dupend <- length(dups)
carbon[dups[dupstart:dupend],'sampleID'] <- "Little Manitwo"
colnames(carbon)[which(colnames(carbon) == "sampleID")] <- "lakeName"
############################################
############################################
############################################

## some data are long format, but with different lengths (nrow) so can't unstack (at least not
##    with my level of R. So, splitting by variableID to be merged....
varsplit <- with(vars, split(vars, variableID))
vars0split <- with(vars0, split(vars0, variableID))
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
allvars <- transform(allvars, windMS = wind * (1000/60/60))
## Lenore pH was a typo - let's change to real value
allvars$pH[allvars$lakeName == "Lenore"] <- 9.3
## Arthur longitude is duplicate of its latitude - let's change to real value
allvars$longitude[allvars$lakeName == "Arthur"] <- -105.443097
## Edouard high pH is an anomaly, Heather checked: "August 2013 was when the pH probe was 
##    malfunctioning and I took additional readings with the handheld meter. 
##    I would use the handheld values, so for Ed, Aug was 9.4."
allvars$pH[allvars$lakeName == "Edouard"] <- 9.4

## sort out df with 0 depths
names <- c(names(vars0split))
varlist <- varnames(vars0split, "dataValue", names, c("variableID", "depth"))
vars0df <- data.frame(lakeName = varlist[[1]]$lakeName)

for (i in seq(varlist)) {
  vars0df <- merge(vars0df, varlist[[i]], by = "lakeName", all = TRUE)
}

vars0df <- vars0df[,-which(colnames(vars0df) == "conductivitySLS")]

## there is a typo for Willow Bunch -- the cond value is actually typed in with the value
##    for salinity, so will calculate the conductivity from the specific conductance
grab <- which(vars0df$lakeName == 'Willow Bunch')
vars0df[grab, 'conductivity'] <- vars0df[grab, "specificConductivity"] * 
  (1 + 0.02 * (vars0df[grab, "temperature"] - 25)) 

## prepare for replace operation and select cols from allvars that match with those in
##    0 depth df
zeronames <- colnames(vars0df)
allsub <- subset(allvars, select = zeronames)

nanumbers <- rowSums(is.na(allsub))
datalost <- allsub[rowSums(is.na(allsub)) >= 2,]
nrow(datalost) # 19; for these sites there is no limno for depth = 1m

## mapply converts numerics to characters and I need to avoid this. Hence take out
##    lakeName since that *is* character and therefore confuses transformations
## (http://stackoverflow.com/questions/15490006/r-as-numeric-matrix)
rownames(allsub) <- allsub$lakeName
allsub$lakeName <- NULL
allsub <- as.matrix(allsub)

rownames(vars0df) <- vars0df$lakeName
vars0df$lakeName <- NULL
vars0df <- as.matrix(vars0df)

## replace NA with the appropriate depth 0 values for cond, pH, salinity, temp:
## Grab names out of the df with the NA values to prep subset for 0 depths
nanames <- datalost$lakeName
replacecols <- which(colnames(allsub) %in% colnames(vars0df))

# order cols to match
reorder <- colnames(vars0df)[order(colnames(vars0df))]
vars0df <- vars0df[,c(reorder)]

replacewith <- which(rownames(vars0df) %in% nanames)
replaceme <- which(rownames(allsub) %in% nanames)

allsub[replaceme,replacecols] <- mapply(replace, allsub[replaceme,replacecols], 
                                         values = vars0df[replacewith, ] )

## put allsub back into allvars
newvals <- which(names(allvars) %in% colnames(allsub))
allvars <- allvars[, -newvals]

allvars <- cbind(allvars, allsub)

## save full dataset for running gasflux on
saveRDS(allvars, "../data/private/heathergasfluxsupp.rds")

