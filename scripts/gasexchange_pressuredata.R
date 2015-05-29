## adapted scripts from bottomoftheheap/...canadian climate data
## we want pressure data, but not avaiable for monthly csv downloads so need also to go hourly for this
## same station as Heather.

## genURLS function

genURLS <- function(id, start, end) {
    years <- seq(start, end, by = 1)
    nyears <- length(years)
    years <- rep(years, each = 12)
    months <- rep(1:12, times = nyears)
    URLS <- paste0("http://climate.weather.gc.ca/climateData/bulkdata_e.html?timeframe=1&Prov=SK&StationID=",
                   id,
                   "&hlyRange=1953-01-30%7C2014-12-31&cmdB1=Go&Year=",
                   years,
                   "&Month=",
                   months,
                   "&Day=27",
                   "&format=csv",
                   "&stationID=",
                   id)
    list(urls = URLS, ids = rep(id, nyears * 12), years = years, months = months)
}

## my test request

datapls <- genURLS(id = 28011, start = 1994, end = 2015) # 264 urls! yay...
head(datapls$urls)

## getData function

getData <- function(stations, folder, verbose = TRUE, delete = TRUE) {
    ## form URLS
    urls <- lapply(seq_len(NROW(stations)),
                   function(i, stations) {
                       genURLS(stations$StationID[i],
                               stations$start[i],
                               stations$end[i])
                   }, stations = stations)
 
    ## check the folder exists and try to create it if not
    if (!file.exists(folder)) {
        warning(paste("Directory:", folder,
                      "doesn't exist. Will create it"))
        fc <- try(dir.create(folder))
        if (inherits(fc, "try-error")) {
            stop("Failed to create directory '", folder,
                 "'. Check path and permissions.", sep = "")
        }
    }
 
    ## Extract the data from the URLs generation
    URLS <- unlist(lapply(urls, '[[', "urls"))
    sites <- unlist(lapply(urls, '[[', "ids"))
    years <- unlist(lapply(urls, '[[', "years"))
    months <- unlist(lapply(urls, '[[', "months"))
 
    ## filenames to use to save the data
    fnames <- paste(sites, years, months, "data.csv", sep = "-")
    fnames <- file.path(folder, fnames)
 
    nfiles <- length(fnames)
 
    ## set up a progress bar if being verbose
    if (isTRUE(verbose)) {
        pb <- txtProgressBar(min = 0, max = nfiles, style = 3)
        on.exit(close(pb))
    }
 
    out <- vector(mode = "list", length = nfiles)
    cnames <- c("Date/Time", "Year", "Month","Day", "Time", "Data Quality",
                "Temp (degC)", "Temp Flag", "Dew Point Temp (degC)",
                "Dew Point Temp Flag", "Rel Hum (%)", "Rel Hum Flag",
                "Wind Dir (10s deg)", "Wind Dir Flag", "Wind Spd (km/h)",
                "Wind Spd Flag", "Visibility (km)", "Visibility Flag",
                "Stn Press (kPa)", "Stn Press Flag", "Hmdx", "Hmdx Flag",
                "Wind Chill", "Wind Chill Flag", "Weather")
 
    for (i in seq_len(nfiles)) {
        curfile <- fnames[i]
 
        ## Have we downloaded the file before?
        if (!file.exists(curfile)) {    # No: download it
            dload <- try(download.file(URLS[i], destfile = curfile, quiet = TRUE))
            if (inherits(dload, "try-error")) { # If problem, store failed URL...
                out[[i]] <- URLS[i]
                if (isTRUE(verbose)) {
                    setTxtProgressBar(pb, value = i) # update progress bar...
                }
                next                             # bail out of current iteration
            }
        }
 
        ## Must have downloaded, try to read file
        ## skip first 16 rows of header stuff
        ## encoding must be latin1 or will fail - may still be problems with character set
        cdata <- try(read.csv(curfile, skip = 16, encoding = "latin1"), silent = TRUE)
 
        ## Did we have a problem reading the data?
        if (inherits(cdata, "try-error")) { # yes handle read problem
            ## try to fix the problem with dodgy characters
            cdata <- readLines(curfile) # read all lines in file
            cdata <- gsub("\x87", "x", cdata) # remove the dodgy symbol for partner data in Data Quality
            cdata <- gsub("\xb0", "deg", cdata) # remove the dodgy degree symbol in column names
            writeLines(cdata, curfile)          # write the data back to the file
            ## try to read the file again, if still an error, bail out
            cdata <- try(read.csv(curfile, skip = 16, encoding = "latin1"), silent = TRUE)
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
 
        ## Must have (eventually) read file OK, add station data
        cdata <- cbind.data.frame(StationID = rep(sites[i], NROW(cdata)),
                                  cdata)
        names(cdata) <- c("StationID", cnames) # ew edited here
        out[[i]] <- cdata
 
        if (isTRUE(verbose)) { # Update the progress bar
            setTxtProgressBar(pb, value = i)
        }
    }
 
    out                                 # return
    saveRDS(out, "outgas.rds") # ew edited here
}

## my station data frame

# nb all stations can cover my year range apart from outlook which has hourly only from 1996... there 
# are randomly monthly data going way back but the monthly interface does not show pressure... if we want to request 
# data it costs.... $100!!!! http://climate.weather.gc.ca/newPriceAnnounce_e.html

#maximum example now
stations <- data.frame(StationID = c(3002, 51441, 3062, 44203, 48977, 3318, 2926, 2925), 
		start = c(1994, 2013, 1994, 2005, 2011, 1996, 1994, 1996), 
		end = c(2013, 2015, 2005, 2011, 2015, 2015, 1996, 2015)) 
			# see weatherstations.csv in /RScripts; these will cover lakes b,c,d,k,l,p,ww.

met <- getData(stations, folder = "./weatherdata")


