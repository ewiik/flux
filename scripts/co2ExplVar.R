## This script creates an rds that contains all relevant data for running regressions against pCO2
## 1. download weather pattern data from web: PDO, SOI, NAO
## 2. process routines sampling data tables

## Part 1:
## ==================================================================================================

## urls, and filenames to use to save the data
urls <- c("http://research.jisao.washington.edu/pdo/PDO.latest", 
          "http://www.cpc.noaa.gov/data/indices/soi")
fileIDs <- c("pdo", "soi")
fnames <- paste0("data/", fileIDs, ".csv")
nfiles <- length(fnames) # I don't think I need file.path here
out <- vector(mode = "list", length = nfiles)

for (i in length(urls)) {
  curfile <- fnames[i]
  # dload <- readLines(urls[1])
  dload <- try(download.file(urls[1], destfile = fnames[1], quiet = TRUE)) 
  cdata <- try(read.csv(fnames[1], skip = 27, encoding = "latin1"), silent = TRUE)
  colnames(cdata) <- cdata[1,] #here it is still data frame
  cdata <- cdata[-c(118:128),] #after this step it becomes a factor...
  cdata <- as.data.frame(cdata) #colnames becomes "cdata", 2015 is not in the year column
  colnames(cdata) <- cdata[1,] #here it is still data frame
  
  
  ## Have we downloaded the file before?
  if (inherits(dload, "try-error")) { # If problem, store failed URL...
    out[[i]] <- URLS[i]
    if (isTRUE(verbose)) {
      setTxtProgressBar(pb, value = i) # update progress bar...
    }
    next                             # bail out of current iteration
  }

  ## Must have downloaded, try to read file
  ## skip first 19 rows of header stuff
  ## encoding must be latin1 or will fail - may still be problems with character set
  cdata <- try(read.csv(curfile, skip = 19, encoding = "latin1"), silent = TRUE)
  cdata <- readLines(fnames[1])
  cdata <- read.table(text = cdata, skip = 28, encoding = "latin1", nrows = nrow(cdata) - 15)
  # !!! THIS IS AS FAR AS I GOT; check http://stackoverflow.com/questions/21192121/
  #                                    r-import-csv-skip-first-and-last-lines
  ## Did we have a problem reading the data?
  if (inherits(cdata, "try-error")) { # yes handle read problem
    ## try to fix the problem with dodgy characters
    cdata <- readLines(curfile) # read all lines in file
    cdata <- gsub("\x87", "x", cdata) # remove the dodgy symbol for partner data in Data Quality
    cdata <- gsub("\xb0", "deg", cdata) # remove the dodgy degree symbol in column names
    writeLines(cdata, curfile)          # write the data back to the file
    ## try to read the file again, if still an error, bail out
    cdata <- try(read.csv(curfile, skip = 19, encoding = "latin1"), silent = TRUE)
      if (inherits(cdata, "try-error")) { # yes, still!, handle read problem
        if (delete) {
        file.remove(curfile) # remove file if a problem & deleting
        }
      out[[i]] <- URLS[i]    # record failed URL...
      if (isTRUE(verbose)) {
        setTxtProgressBar(pb, value = i) # update progress bar...
        }
    next                  # bail out of current iteration
      }
  }

  ## Must have (eventually) read file OK
  out[[i]] <- cdata

  if (isTRUE(verbose)) { # Update the progress bar
  setTxtProgressBar(pb, value = i)
    }
  out  # return
  }

## Part 2:
## ==================================================================================================

## read in supporting monitoring data tables grabbed from database
## I created a Date2 column in OpenOffice to convert the current display of month as character to month
##    as numeric (former is how it came from database into excel cause excel sucks)

routines <- read.csv("data/private/qpco2supportdata.csv")
routines <- transform(routines, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))
routines <- routines[with(routines, order(LAKE, Date)),]

produc <- read.csv("data/private/Production.csv")
produc <- transform(produc, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))

chl <- read.csv("data/private/Chl.csv")
chl <- transform(chl, Date = as.POSIXct(as.character(Date2), format = "%Y-%m-%d"))

## remove extra date column (originally retained in case wanna check that the date format
##    conversion worked in OpenOffice)
routines <- routines[,-which(colnames(routines)=="Date2")] # wanted to avoid using numbers for columns
produc <- produc[,-which(colnames(produc)=="Date2")] 
chl <- chl[,-which(colnames(chl)=="Date2")] 


## split produc and chl for getting means for replicated production estimates
prodsplit <- with(produc, split(produc, list(LAKE, Date), drop = TRUE))
chlsplit <- with(chl, split(chl, list(LAKE, Date), drop = TRUE))

## create wrapper that does colmeans for the columns we want and spits out df we can merge
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(LAKE = df['LAKE'][1,], Date = df['Date'][1,]), t(means))
}

## choose columns we want means for
## FIXME: confirm with kerri that we want netoxy, and which respiration measure we want
## FIXME: confirm with kerri whether we want integrated or surface chlorophyll!
## FIXME: chl columns: if we take "total chl", there are NAs where method only measured chla ...
##    so need to know if we are taking total or chl a ... and if we can replace missing data

prodmeans <- lapply(prodsplit, mycolMeans, cols = c("NetOxy_ppm", "Resp_mgC_m3_H")) 
prodmeans <- do.call(rbind, prodmeans)
rownames(prodmeans) <- NULL
prodmeans <- prodmeans[with(prodmeans, order(LAKE, Date)),]

chlmeans <- lapply(chlsplit, mycolMeans, cols = c("CHLC", "Total_chl")) # once chl issue is fixed, we 
#   will need to change this and subset by TreatmentNewLabel (or Old... check if no difference)
chlmeans <- do.call(rbind, chlmeans)
rownames(chlmeans) <- NULL
chlmeans <- chlmeans[with(chlmeans, order(LAKE, Date)),]

## merge all dfs by LAKE and Date
co2explained <- merge(chlmeans, prodmeans) ## FIXME: NAs are now NaN.. problem for later?
co2explained <- merge(co2explained, routines)

## save output for later
saveRDS(co2explained, "data/private/co2explained.rds")
