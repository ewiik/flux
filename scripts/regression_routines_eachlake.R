## following from regression_routines -pH and -models, this script creates lake-
##    specific models following meeting with Peter 25.01.2015
## final models chosen were co2mod and its residuals as resmodred for CO2 ~.
##    and egmodred.2 for pH ~ . These models are in this script run, the others
##    are set to not run

## set to whether or not to run the discarded models, and whether or not to plot output
##    from models that were kept
runextras <- FALSE
plotmods <- FALSE

## load necessary packages
library("mgcv")
library("ggplot2")

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
if(runextras) {
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

## save summary as txt document
egmodsum <- summary(egmod)
sink("../data/private/egmodsummary.txt")
egmodsum
sink()
}
## 2. In fact, good to log some of the more skewed distributions
## change -ve values to 0, then add 1
regvarf2 <- regvarf
regvarf2$`Chl_a_ug_L`[regvarf2$`Chl_a_ug_L` <=0 ] <- 0
regvarf2$`DOC_mg_L`[regvarf2$`DOC_mg_L` <=0 ] <- 0
regvarf2$`Chl_a_ug_L` <- regvarf2$`Chl_a_ug_L` + 1
regvarf2$`DOC_mg_L` <- regvarf2$`DOC_mg_L` + 1
regvarf2 <- transform(regvarf2, dummy = rep(1, nrow(regvarf2)))

## for this model, we are removing some of the Lake-specificity where nothing signif
##    was found in previous model
if(runextras) {
egmod.red <- gam(pH_surface ~ 
                   s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                   s(GPP_h) + s(GPP_h, by = Lake, m = 1) +
                   s(log10(TDN_ug_L)) + 
                   s(log10(DOC_mg_L)) +
                   s(Oxygen_ppm) + s(Oxygen_ppm, by = Lake, m = 1) +
                   te(PDO, SOI) +
                   s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator),
                 data = regvarf2,
                 select = TRUE, method = "REML", family = gaussian,
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE))
                 # this setting makes the computation a bit more efficient 
## longer object length is not a multiple of shorter object length!!??

egmodredsum <- summary(egmod.red)
sink("../data/private/egmodredsummary.txt")
egmodredsum
sink()
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
                    s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator
                  data = regvarf2,
                  select = TRUE, method = "REML", family = scat(),
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, trace = TRUE,
                                        newton = list(maxHalf = 60)))
egmodred2sum <- summary(egmod.red2)
if (plotmods) {
sink("../data/private/egmodred2summary.txt")
egmodred2sum
sink()
}

## plot output of model
if (plotmods) {
pdf("../docs/egmod-reduced-log-covariates.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(egmod.red2, pages = 4, scheme = 2)
par(op)
dev.off()
}

## Get site-specific data/predictions
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

chla.pdat <- with(droplevels(regvarf2),
                  data.frame(`Chl_a_ug_L` = rep(seq(min(`Chl_a_ug_L`, na.rm = TRUE),
                                                    max(`Chl_a_ug_L`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
chla.pdat <- merge(chla.pdat, lakeXbar)
chla.pred <- predict(egmod.red2, newdata = chla.pdat, type = "terms")
whichCols <- grep("Chl", colnames(chla.pred))
chla.pdat <- cbind(chla.pdat, Fitted = rowSums(chla.pred[, whichCols]))

ggplot(chla.pdat, aes(x = Chl_a_ug_L, y = Fitted, colour = 
                        ifelse(Lake == "WW", "Wascana", "others"))) +
  geom_line() + 
  scale_colour_manual(name = 'Lakes', values = c("#999999", "#E69F00")) +
  xlab('Chl a (ug/L)') + ylab('mean-0 pH')
## see http://stackoverflow.com/questions/11838278/
##    plot-with-conditional-colors-based-on-values-in-r

  
## use final models created for CO2, adapted to see if any one Lake deviates
##    from the global trend
## =================================================================================
## 1. can we use a blanket model for CO2 based on pH?
co2mod <- gam(co2Flux ~ 
               s(pH_surface) + s(pH_surface, by = Lake, m = 1) +
               # m = 1 means that penalty goes with 1st derivative... required
               #    here because we're doing both s() and s(..,by = ...)
               s(Lake, Year, bs = "re"), # still need to explicitly keep Lake as re!!!
             data = regvarf,
             select = TRUE, method = "REML", family = gaussian,
             na.action = na.exclude,
             control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))

co2modsum <- summary(co2mod)
if (plotmods) {
sink("../data/private/co2modsummary.txt")
co2modsum
sink()
}

gam.check(co2mod)
## FIXME: the kindex and pvalues are still NA for this model too

if(runextras) {
co2modnull <- gam(co2Flux ~ 
                s(pH_surface) +
                s(Lake, Year, bs = "re"), 
              data = regvarf,
              select = TRUE, method = "REML", family = gaussian,
              na.action = na.exclude,
              control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))

anova(co2modnull, co2mod, test = "LRT")
AIC(co2modnull, co2mod)
# --> we need to take the more complex model
}

## plot output of model
if (plotmods) {
pdf("../docs/co2mod.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(co2mod, pages = 4, scheme = 2)
par(op)
dev.off()
}

## model residuals
res <- resid(co2mod, type = "pearson")

if(runextras) {
resmod <- gam(res ~ 
                    s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                    s(GPP_h) + s(GPP_h, by = Lake, m = 1) +
                    s(log10(TDN_ug_L)) + s(log10(TDN_ug_L), by = Lake, m = 1) +
                    s(log10(DOC_mg_L)) + s(log10(DOC_mg_L), by = Lake, m = 1) +
                    s(Oxygen_ppm) + s(Oxygen_ppm, by = Lake, m = 1) +
                    te(PDO, SOI) +
                    s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator),
                  data = regvarf2,
                  select = TRUE, method = "REML", family = scat(),
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, trace = TRUE,
                                        newton = list(maxHalf = 60)))

co2resmodsum <- summary(resmod)
sink("../data/private/co2resmodsummary.txt")
co2resmodsum
sink()

pdf("../docs/resmod-full-covariates.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(resmod, pages = 5, scheme = 2)
par(op)
dev.off()
}

## some of the GPP lake-specific things look a bit far-fetched so took Lake away there.
resmodred <- gam(res ~ 
                   s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                   s(GPP_h) +
                   s(log10(TDN_ug_L)) + 
                   s(log10(DOC_mg_L)) + 
                   s(Oxygen_ppm) + 
                   te(PDO, SOI) +
                   s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator
                 data = regvarf2,
                 select = TRUE, method = "REML", family = scat(),
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE,
                                       newton = list(maxHalf = 60)))
if (runextras) {
anova(resmodred, resmod, test = "LRT")
## the simpler model is just fine
}

co2resmodredsum <- summary(resmodred)

if (plotmods) {
sink("../data/private/co2resmodredsummary.txt")
co2resmodredsum
## s(log10(TDN_ug_L))          2.931e+00      9 25.756 9.97e-06 ***
sink()

plot(resmodred, pages = 3, scheme = 2)
}
## can we simplify even further?
if (runextras) {
resmodredchl <- gam(res ~ 
                      s(log10(Chl_a_ug_L)) +
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

anova(resmodredchl, resmodred, test = "LRT")
## better to keep the more complicated model
}

### Get site specific data/predictions
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

chla.pdat <- with(droplevels(regvarf2),
                  data.frame(`Chl_a_ug_L` = rep(seq(min(`Chl_a_ug_L`, na.rm = TRUE),
                                                    max(`Chl_a_ug_L`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
chla.pdat <- merge(chla.pdat, lakeXbar)
chla.pred <- predict(resmodred, newdata = chla.pdat, type = "terms")
whichCols <- grep("Chl", colnames(chla.pred))
chla.pdat <- cbind(chla.pdat, Fitted = rowSums(chla.pred[, whichCols]))

ggplot(chla.pdat, aes(x = Chl_a_ug_L, y = Fitted, colour = Lake)) +
  geom_line()



if (runextras) {
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
}
