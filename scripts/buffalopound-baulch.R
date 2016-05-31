### quick look at Buffalo Pound monitoring data from 1985
## base file based on csv in frompeter/

## read files and subset
if(!file.exists("../data/private/MasterfileBuffaloPoundWaterQualityData15032016.csv")) {
  stop('get BUffalo Pound weekly survey data from Peter or Emma')
}
buff <- read.csv("../data/private/MasterfileBuffaloPoundWaterQualityData15032016.csv")
remstart <- which(names(buff)=='Anacystis')
buffsub <- buff[,-c(remstart:length(names(buff)))]

if(!file.exists('../data/maunaloa.csv')) {
  source('../functions/getmaunaloa.R') }
ml <- read.csv('../data/maunaloa.csv')

## read in salcalc and co2 flux function
source('../data/private/salcalc.R')
source("../functions/gasExchangeExtra.R")

## calculate salinity from temperature and other things
## remember function converts from uS to mS so we don't need to change this
## assuming constant pressure here as for pham data: i.e. 95kPa
buffsub$dbar <- 95/10 + 1
buffsub$salcalc <- with(buffsub, salcalc(temp = temp, cond = cond, dbar=dbar))
saltozero <- which(buffsub$salcalc < 0) # some -ves
buffsub[saltozero, 'salcalc'] <- 0

## some 0 dic should be na (bicarb and carb as na)
dicna <- which(buffsub$bicarb == 0)
buffsub[dicna,c('bicarb','carb')] <- NA

## do TIC as per pham's stuff (see pham-DICcalcs.R)
## proportion of (H)CO3 weight that is C * total weight
## molecular weights of ions as g/mol:
bicarbweight <- 1.0079 + 12.0107 + 3*15.9994 # g/mol 
carbweight <- 12.0107 + 3*15.9994 # g/mol 
co2weight <- 12.0107 + 2*15.9994 # g/mol
# use molecular weights to calculate mg/L C
buffsub$bicarbC <- (12.0107/bicarbweight) * buffsub$bicarb
buffsub$carbC <- (12.0107/carbweight) * buffsub$carb

## calculate k constants (from Kerri's gas flux calculations)
buffsub$ionstrength <- with(buffsub, cond*7*0.0000025)

buffsub$pk1 <- with(buffsub, (0.000011*temp^2 - 0.012*temp + 6.58) - (0.5*sqrt(ionstrength))/
                    (1 + 1.4*sqrt(ionstrength)))
buffsub$pk2 <- with(buffsub, (0.00009*temp^2 - 0.0137*temp + 10.62) - (2*sqrt(ionstrength))/
                    (1 + 1.4*sqrt(ionstrength)))

buffsub$k1 <- with(buffsub, 10^(-pk1))
buffsub$k2 <- with(buffsub, 10^(-pk2))

## calculate H+ concentration from measured pH
buffsub$hplus <- with(buffsub, 10^(-ph))

## calculate a0-2 as per Kerri's gas flux calculations for 
##    proportions of DIC as CO2, HCO3 and CO3 respectively
buffsub$ao <- with(buffsub, 1/(1 + 10^(-pk1 + ph) + 10^(-(pk1 + pk2) + 2*ph)))
buffsub$a1 <- with(buffsub, 1/(10^(-ph+pk1) + 1 + 10^(-pk2 + ph)))
buffsub$a2 <- with(buffsub, 1/(10^(-2*ph + (pk1 + pk2)) + 10^(-ph + pk2) + 1))

## null way just to check
buffsub$totaldicnull <- with(buffsub, bicarbC + carbC)
## then pham way
buffsub$totaldic <- with(buffsub, (bicarbC+carbC)/(a1+a2))
## one outlier due to a pH glitch in the data
check <- which(buffsub$totaldic > 5000)
buffsub[check,'ph'] <- NA
buffsub[check,'totaldic'] <- buffsub[check,'totaldicnull']

## found another outlier
phtona <- which(buffsub$ph > 15) #it's 25!!
buffsub[phtona, 'ph'] <- NA
## change units that need changing
buffsub$dicumol <- buffsub$totaldic/0.012 

## and for cond...
condtona <- which(buffsub$cond < 10)
buffsub[condtona, 'cond'] <- NA

## merge with mauna loa data
buffsub <- merge(buffsub, ml, by.x = c('year', 'month'), by.y = c('Year','Month'))
names(buffsub)[which(names(buffsub) == 'pCO2')] <- 'pCO2atm'

## calculate flux
## FIXME: what wind to plug in? pco2atm? now 2004-2004 approximate early summer value
gasfluxes <- with(buffsub, gasExchangeExtra(temp = temp, cond = cond, ph=ph, wind=4, 
                                         alknotdic = FALSE, dic = dicumol, 
                                         salt = salcalc, kpa=95, pco2atm = pCO2atm))
with(gasfluxes,plot(density(pco2, na.rm=TRUE)))

buffsub <- cbind(buffsub, gasfluxes)

## plot stuff
with(buffsub, plot(pco2 ~ ph))
abline(400,0)
with(buffsub, plot(fluxenh ~ ph))
abline(0,0)

## rename to make data sheet
names(buffsub)[which(names(buffsub) %in% c('fluxenh', 'pco2'))] <-
  c("CO2fluxmmolm2d", "pCO2uatm")

## save
saveRDS(buffsub, '../data/private/bp-longterm-mod.rds')
