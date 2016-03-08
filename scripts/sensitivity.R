## sensitivity analysis as to how much pH dominates the CO2 flux of our study sites
## see https://cran.r-project.org/web/packages/pse/vignettes/pse_tutorial.pdf
##    https://cran.r-project.org/web/packages/sensitivity/sensitivity.pdf
## what indices to be used for expressing sensitivity?

## source necessary packages:
library("mc2d")
library("pse")
library("fitdistrplus")

## source the function we want to decompose
source("../functions/gasExchangeFlex.R")

## source files with parameters
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")

## define the parameters
The 1st thing that must be done is to determine exactly what are the input pa-
  rameters to your model.  You should list which parameters should be investigated,
what are the probability density functions (PDFs) from which the parameter values
will be calculated, and what are the arguments to these PDFs.
factors <- c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", "TICumol", "Pressure")
distro <- c("qnorm", "qnorm", "qunif")
props <- list( list(mean=1.7, sd=0.3), list(mean=40, sd=1),
                 +         list(min=1, max=50) )

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
# non-normal: salcalc, meanwind, pressure, TIC, pH, cond, temp
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

## see list http://www.statmethods.net/advgraphs/probability.html
## decide which to look at and what likely dist to test with
## FIXME: not working very well for t distribution! 
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
x <- fluxes$SalCalc
disttest <- "weibull"
## chisq may require this to run: start = list(df = 8) or some other number
## t retuires same to run but behaves suspiciously (see plots with expected lines)
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
