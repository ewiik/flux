## following from regression_routines -pH and -models, this script creates lake-
##    specific models following meeting with Peter 25.01.2015

## load necessary packages
library("mgcv")

## load necessary data
if (!file.exists("../data/private/regvars.rds")) {
  source("../scripts/regression_routines.R")
}
regvars <- readRDS("../data/private/regvars.rds")
regvarf <- regvars
regvarf <- transform(regvarf, Year = as.factor(Year)) # make Year into factor for re

## split regvarf by lake then remove all NAs from list
regsplit <- with(regvarf, split(regvarf, list(Lake)))
regsplit <- regsplit[lapply(regsplit,nrow) > 1]

## use final models created for CO2, adapted to see if any one Lake deviates
##    from the global trend
## =================================================================================
## 1. model most similar to ones done before
egmod <- gam(pH_surface ~ 
               s(Chl_a_ug_L) + s(Chl_a_ug_L, by = Lake, m = 1) +
               # m = 1 means that penalty goes with 1st derivative... required
               #    here because we're doing both s() and s(..,by = ...)
               s(GPP_h) + s(GPP_h, by = Lake, m = 1) +
               s(TDN_ug_L) + s(TDN_ug_L, by = Lake, m = 1) + 
               s(DOC_mg_L) + s(DOC_mg_L, by = Lake, m = 1) +
               s(Oxygen_ppm) + s(Oxygen_ppm, by = Lake, m = 1) +
               te(PDO, SOI) + te(PDO, SOI, by = Lake, m = 1) +
               s(Lake, Year, bs = "re"), # still need to explicitly keep Lake as re!!!
             data = regvarf,
             select = FALSE, method = "REML", family = gaussian,
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE))

summary <- FALSE
if(summary) {
  summary(egmod)
}

## 2. In fact, good to log some of the more skewed distributions
## change -ve values to 0, then add 1
regvarf2 <- regvarf
regvarf2$`Chl_a_ug_L`[regvarf2$`Chl_a_ug_L` <=0 ] <- 0
regvarf2$`DOC_mg_L`[regvarf2$`DOC_mg_L` <=0 ] <- 0
regvarf2$`Chl_a_ug_L` <- regvarf2$`Chl_a_ug_L` + 1
regvarf2$`DOC_mg_L` <- regvarf2$`DOC_mg_L` + 1

## for this model, we are removing some of the Lake-specificity where nothing signif
##    was found in previous model
egmod.red <- gam(pH_surface ~ 
                   s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                   s(GPP_h) + s(GPP_h, by = Lake, m = 1) +
                   s(log10(TDN_ug_L)) + 
                   s(log10(DOC_mg_L)) +
                   s(Oxygen_ppm) + s(Oxygen_ppm, by = Lake, m = 1) +
                   te(PDO, SOI) +
                   s(Lake, Year, bs = "re"),
                 data = regvarf2,
                 select = TRUE, method = "REML", family = gaussian,
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE))
                 # this setting makes the computation a bit more efficient 

summary <- FALSE
if(summary) {
  summary(egmod.red)
}

## 3. model that removes further unnecessary elements to make fitting faster and simpler
##    AND introduces scat() again. Testing to see if the step failure can be averted
##    by setting maxHalf = 60
egmod.red2 <- gam(pH_surface ~ 
                    s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                    s(GPP_h) +
                    s(log10(TDN_ug_L)) + 
                    s(log10(DOC_mg_L)) +
                    s(Oxygen_ppm) +
                    te(PDO, SOI) +
                    s(Lake, Year, bs = "re"),
                  data = regvarf2,
                  select = TRUE, method = "REML", family = scat(),
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, trace = TRUE,
                                        newton = list(maxHalf = 60)))
summary <- TRUE
if(summary) {
  summary(egmod.red2)
}

## plot output of model
pdf("../docs/egmod-reduced-log-covariates.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(egmod.red2, pages = 4, scheme = 2)
par(op)
dev.off()

## ===========================================================================
## Below are my initial apply applications before learning that I can introduce
##    by = ... into a gam call which makes it all better cause I am still using 
##    all the data that way and other obvious reasons
## ===========================================================================
co2gams <- function(df) {
  if (nrow(df) < 3) { return(NA) } else {
    co2gam <- gam(co2Flux ~ s(pH_surface, k = 20), data = df,
             select = TRUE, method = "ML", family = scat(),
             na.action = na.exclude)
  }}

resgams <- function(df, res) {
  if (nrow(df) < 3) { return(NA) } else {
  res <- as.vector(res)
    resmod <- gam(res ~ s(Chl_a_ug_L) + s(GPP_h) + s(TDN_ug_L) + 
             s(DOC_mg_L) + s(Oxygen_ppm) + 
             te(PDO, SOI) + s(Year, bs = "re"), data = df,
           select = TRUE, method = "REML", family = scat(),
           na.action = na.exclude) }
  }

## apply gam 
regmods <- lapply(regsplit, co2gams)

## apply residual extraction to regmods
allres <- lapply(regmods, resid, type = "pearson", na.action = na.exclude)

## apply resgams to the residuals
## FIXME: this not doing the right thing!
resmods <- mapply(resgams, regsplit, allres)

## can lapply these to the above to check diagnostics however
##    over gam.check the plotting default means I don't know how to store objects
lapply(regmods, gam.check)
regsumms <- lapply(regmods, summary)
## deviance explained ranges from 68 - 80%
lapply(regsumms, '[[', "dev.expl")

lapply(regmods, plot, pers = TRUE, pages = 1)


## 2. pH ~ .
## ===========================================================================
phgams <- function(df) {
  if (nrow(df) < 3) { return(NA) } else {
  gam(pH_surface ~ s(Lake, Year, bs = "re") + s(Chl_a_ug_L) + s(GPP_h) +
               s(TDN_ug_L) + s(DOC_mg_L) + s(Oxygen_ppm) + te(PDO, SOI), data = df,
             select = TRUE, method = "REML", family = gaussian, na.action = na.exclude)
}}

phmods <- lapply(regsplit, phgams)

summary(phmod)
plot(phmod, pages = 1, pers = TRUE)
gam.check(phmod)

