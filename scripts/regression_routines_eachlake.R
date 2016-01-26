## following from regression_routines -pH and -models, this script creates lake-
##    specific models following meeting with Peter 25.01.2015

## load necessary packages
library("mgcv")

## load necessary data
if (!file.exists("../data/private/regvars.rds")) {
  source("../scripts/regression_routines.R")
}
regvars <- readRDS("../data/private/regvars.rds")

## use final models created for CO2, though lapply this by lake (and therefore
##    remove Lake as a random error)
co2gams <- function(df) {
  if (nrow(df) < 3) { return(NA) } else {
    co2gam <- gam(co2Flux ~ s(pH_surface, k = 20), data = df,
             select = TRUE, method = "ML", family = scat(),
             na.action = na.exclude)
  }}

resgams <- function(df, res) {
  if (nrow(df) < 3) { return(NA) } else {
  res <- as.data.frame(res)
    resmod <- gam(res ~ s(Chl_a_ug_L) + s(GPP_h) + s(TDN_ug_L) + 
             s(DOC_mg_L) + s(Oxygen_ppm) + 
             te(PDO, SOI) + s(Year, bs = "re"), data = df,
           select = TRUE, method = "REML", family = scat(),
           na.action = na.exclude) }
  }

## split regvars by lake then remove all NAs from list
regsplit <- with(regvars, split(regvars, list(Lake)))
regsplit <- regsplit[lapply(regsplit,nrow) > 1]

## apply gam 
regmods <- lapply(regsplit, co2gams)

## apply residual extraction to regmods
allres <- lapply(regmods, resid, type = "pearson", na.action = na.exclude)

## apply resgams to the residuals
resmods <- mapply(resgams, regsplit, allres)
## FIXME: due to this error objects aren't created:
##    Error in diag(c(rep(S11, rank), rep(0, p - rank))) : vector is too large 
## with later model where res is made into data frame, I get this:
##    Error in model.frame.default(formula = res ~ 1 + Chl_a_ug_L + GPP_h +  : 
##    invalid type (list) for variable 'res' 

## can lapply these to the above to check diagnostics however
##    over gam.check the plotting default means I don't know how to store objects
lapply(regmods, gam.check)
regsumms <- lapply(regmods, summary)
## deviance explained ranges from 68 - 80%
lapply(regsumms, '[[', "dev.expl")

lapply(regmods, plot, pers = TRUE, pages = 1)
