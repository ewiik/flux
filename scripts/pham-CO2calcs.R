## this script calculates CO2 for sam's data set

## read in DIC anc cond-complete pham data set
## note that wind, pressure and temp are just plugged in, not determined!
pham <- read.csv('../data/private/phamdat-mod.csv')

## read in CO2 flux calc
source("../functions/gasExchangeExtra.R")

## change units that need changing
pham$dicumol <- pham$totaldic3/0.012 # choosing the added HCO3, CO3 (measured) / constants
#   method with co2 added

## calculate flux
## FIXME: check local mods to gasExchangeFlex before pushing to git
## FIXME: what wind to plug in? pco2atm? now 2004-2004 approximate early summer value
gasfluxes <- with(pham, gasExchangeExtra(temp = temp, cond = condUS, ph=pH_ave, wind=4, 
                                        alknotdic = FALSE, dic = dicumol, 
                                        salt = salsurf_ppt, kpa=95, pco2atm = 380))
with(gasfluxes,plot(density(pco2, na.rm=TRUE)))

pham <- cbind(pham, gasfluxes)
with(pham, plot(pco2 ~ pH_ave, cex = totaldic3/100))
legend('topright', legend = c('DIC = 100mg/L', 'DIC = 400mg/L'), pch=c(1,1), pt.cex=c(1,4))

## rename some column names to have units:
names(pham)[which(names(pham) %in% c('fluxenh', 'pco2'))] <-
  c("CO2fluxmmolm2d", "pCO2uatm")

## save for later
write.csv(pham, '../data/private/phamdat-co2.csv')
