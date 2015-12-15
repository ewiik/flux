## sensitivity analysis as to how much pH dominates the CO2 flux of our study sites
## see https://cran.r-project.org/web/packages/pse/vignettes/pse_tutorial.pdf
##    https://cran.r-project.org/web/packages/sensitivity/sensitivity.pdf

## source necessary packages:
library("mc2d")
library("pse")

## source the function we want to decompose
source("functions/gasExchangeFlex.R")

## source files with parameters
regvars <- readRDS("data/private/regvars.rds")
fluxes <- readRDS("data/private/params-flux.rds")

## define the parameters
## see this: http://stats.stackexchange.com/questions/16646/
##    what-is-a-good-index-of-the-degree-of-violation-of-normality-and-what-descriptive
The rst thing that must be done is to determine exactly what are the input pa-
  rameters to your model.  You should list which parameters should be investigated,
what are the probability density functions (PDFs) from which the parameter values
will be calculated, and what are the arguments to these PDFs.
factors <- c("Temperature", "Conductivity", "pH", "meanWind", "SalCalc", "TIC", "Pressure")
distro <- c("qnorm", "qnorm", "qunif")
props <- list( list(mean=1.7, sd=0.3), list(mean=40, sd=1),
                 +         list(min=1, max=50) )

qqplott <- function(vect) {
  qqnorm(vect, ylab = mean(vect, na.rm = TRUE))
  qqline(vect)
}
apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWind", "SalCalc", "TIC", "Pressure")], 
      2, qqplott)   
apply(fluxes[,c(7:8,11,13,17:18,22)], 2, shapiro.test)   
# normal: 
# non-normal: salcalc, meanwind, pressure, TIC, pH, cond, temp