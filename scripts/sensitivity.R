## sensitivity analysis as to how much pH dominates the CO2 flux of our study sites
## see https://cran.r-project.org/web/packages/pse/vignettes/pse_tutorial.pdf
##    https://cran.r-project.org/web/packages/sensitivity/sensitivity.pdf
## what indices to be used for expressing sensitivity?
## resources:
##  1. visual distribution cheat sheet http://www.itl.nist.gov/div898/handbook/eda/section3/eda366.htm
##  2. list of R shorts for distributions http://www.statmethods.net/advgraphs/probability.html
##  3. blurb on chisq: http://www.civil.uwaterloo.ca/brodland/EasyStats/EasyStats/Chi_squared_Distribution.html
##  4. fitdist issues: http://stats.stackexchange.com/questions/158163/why-does-this-data-throw-an-error-in-r-fitdistr
##  5. fitting distros: http://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
##  6. chisq: https://www.youtube.com/watch?v=hcDb12fsbBU
##          http://www.civil.uwaterloo.ca/brodland/EasyStats/EasyStats/Chi_squared_Distribution.html


## source necessary packages:
library("mc2d")
library("pse")
library("fitdistrplus") # loads MASS too

## source the function we want to decompose
source("../functions/gasExchangeFlex.R")

## source files with parameters
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")

## define the parameters
factors <- c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", "TICumol", "Pressure")

## Get distributions for each variables
## ====================================
## basic distro plotting
dplott <- function(vect) {
  plot(density(vect, na.rm = TRUE), ylab = mean(vect, na.rm = TRUE))
}
opar <- par(mfrow = c(2,4))
apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure")], 2, dplott) 
par(opar)

qqplott <- function(vect) {
  qqnorm(vect, ylab = mean(vect, na.rm = TRUE))
  qqline(vect)
}
opar <- par(mfrow = c(2,4))
apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure")], 
      2, qqplott)  
par(opar)
apply(fluxes[,c(7:8,11,13,17:18,22)], 2, shapiro.test)   
# normal: NONE
## see this: http://stats.stackexchange.com/questions/16646/
##    what-is-a-good-index-of-the-degree-of-violation-of-normality-and-what-descriptive

## So let's try the packages on fitting distributions
giveDist <- function(vect) {
  if(any(is.na(vect))) {
    remove <- which(is.na(vect))
    descdist(vect[-remove], discrete = FALSE) }
  else {descdist(vect, discrete = FALSE)}
  legend("bottomleft", legend = head(vect))
}
opar <- par(mfrow = c(2,4))
apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure")], 
      2, giveDist) 
par(opar)
## --> normal:  dicumol, salcalc, (meanwind)
## --> unif: cond, pressure
## normal-logistic: pH               
## bleurgh: temperature

## Look at uncertain-distro variables and test with likely distro 
## "For the following named distributions, reasonable starting values will be computed 
##    if start is omitted (i.e. NULL) : "norm", "lnorm", "exp" and "pois", "cauchy", 
##    "gamma", "logis", "nbinom" (parametrized by mu and size), "geom", "beta", 
##    "weibull" from the stats package; "invgamma", "llogis", "invweibull", "pareto1", 
##    "pareto" from the actuar package. Note that these starting values may not be 
##    good enough if the fit is poor. The function uses a closed-form formula to fit the 
##    uniform distribution. If start is a list, then it should be a named list with the 
##    same names as in the d,p,q,r functions of the chosen distribution. 
##    If start is a function of data, then the function should return a named list with 
##    the same names as in the d,p,q,r functions of the chosen distribution. 
x <- fluxes$Pressure
disttest <- "weibull"
## chisq may require this to run: start = list(df = 8) or some other number
## t requires same to run but behaves suspiciously (see plots with expected lines)
{if(any(is.na(x))) {
  remove <- which(is.na(x)) 
  fit.test <- fitdist(x[-remove], disttest)
  fit.norm <- fitdist(x[-remove], "norm")
}
else { fit.test <- fitdist(x, disttest)
fit.norm <- fitdist(x, "norm")}
}

plot(fit.test)
plot(fit.norm)

## lower = better, not in the absolute value sense
fit.test$aic
fit.norm$aic
## for temperature, weibull > normal > chisq | gamma --> weibull
## for pH, logis > normal > cauchy
## for cond, pressure: unif doesn't produce AIC but by eye, normal is better
## for DICumol, normal > logis > cauchy
## looking at dist plot, normal for windspeed
## for SalCalc, weibull > normal

## Continue defining parameters for the chosen distros
distro <- c("qweibull", "qnorm", "qlogis", "qnorm", "qweibull", "qnorm", "qnorm") 
# in order of 'factors'; retrieve args to list by calling the fit.[] object
props <- list( list(shape=4.02, scale=18.65), list(mean=1013, sd=499), 
               list(location=8.84,scale=0.31), list(mean=4.98, sd=0.83),
               list(shape=2.13, scale=0.68), list(mean=3821, sd=1068),
               list(mean=94.6, sd=0.17))

## Consider interrelationships between parameters
## FIXME: define this numerically once find out the best way to achieve
pairs(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure")])
## WHAT T ~ meanWind!! (though a lot of noise)
## cond ~ SalCalc ~ TICumol --> these need a rank order combination

# "The  model  that  you  wish  to  analyse  must  be  formulated  as  an R function  that
# receives a data.frame in which every column represent a different parameter, and
# every line represents a different combination of values for those parameters.  The
# function must return an array with the same number of elements as there were lines
# in the original data frame,  and each entry in the array should correspond to the
# result  of  running  the  model  with  the  corresponding  parameter  combination"
## CHECK!