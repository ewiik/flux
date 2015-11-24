## use gasexchange_precipdata.R to source regina info for routines analysis and
##    make info suitable to combine with other variables

## source necessary functions:
source("functions/gasexchange_precipdata.R")

## my station data frame: for regina gilmour as per heather's searches
stations <- data.frame(StationID = c(3007), start = c(1994), end = c(2014))

## create precipitation data frame
precip <- getData(stations, folder = "data/precipdata/routines")
colvector <- c("StationID", "Date/Time", "Year", "Month", 'Day', "Total Rain (mm)", "Total Snow (cm)", 
               "Total Precip (mm)", "Snow on Grnd (cm)")
precips <- lapply(precip, "[", colvector) # need to circumvent fact that rbind won't roll on odd column 
#   contents that result in conflicting column class assignments for different dfs in the list

precips <- do.call("rbind", precips)
rownames(precips) <- NULL

saveRDS(precips, "data/precip.rds") # save for future

## hmmm... annual snow total? wrapper for which col want
mycolSums <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  sum <- colSums(subdf, na.rm = TRUE)
  cbind(data.frame(Year = df['Year'][1,], t(sum)))
}

precipsplit <- with(precips, split(precips, list(Year), drop = TRUE))
snowtotal <- lapply(precipsplit, mycolSums, cols = c("Total Snow (cm)")) 
snowtotal <- do.call("rbind", snowtotal)
rownames(snowtotal) <- NULL
