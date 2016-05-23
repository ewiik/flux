## calculating DIC from other variables
## we've got pH, hco3 and co3.Note that HCO3 and CO3 are in mg/L total ion, not C

## load packages
library('oce')

## read in file and add variables that we need
## FIXME: was temperature available somewhere?
## Kerri - we decided to use 20 deg to be able to compare across lakes 
testing <- read.csv('../data/private/phamdat.csv')
testing$temp <- rep(20) # assuming 20
testing$pressure <- rep(1) # assuming 1atm since we are at water surface

## remove NAs
remove <- which(is.na(testing$salsurf_ppt))
nonas <- testing[-remove,]

## remove Wolverine which has the wrong pH in
makena <- which(nonas$pH_ave < 6)
nonas$pH_ave[makena] <- NA

## convert Pham's ion data to mg/L in Carbon --> 
##    proportion of (H)CO3 weight that is C * total weight
## molecular weights of ions as g/mol:
bicarbweight <- 1.0079 + 12.0107 + 3*15.9994 # g/mol 
carbweight <- 12.0107 + 3*15.9994 # g/mol 
co2weight <- 12.0107 + 2*15.9994 # g/mol
# use molecular weights to calculate mg/L C
nonas$bicarbC <- (12.0107/bicarbweight) * nonas$Bicarb_mgL
nonas$carbC <- (12.0107/carbweight) * nonas$carb_mgL

## calculate conductivity from salinity
nonas$rratio <- with(nonas, swCSTp(salsurf_ppt, temp, pressure))
nonas$condMS <- with(nonas, rratio*42.914) # milliSiemens
nonas$condUS <- with(nonas, condMS * 1000) # microSiemens

## calculate k constants (from Kerri's gas flux calculations)
nonas$ionstrength <- with(nonas, condUS*7*0.0000025)

nonas$pk1 <- with(nonas, (0.000011*temp^2 - 0.012*temp + 6.58) - (0.5*sqrt(ionstrength))/
                    (1 + 1.4*sqrt(ionstrength)))
nonas$pk2 <- with(nonas, (0.00009*temp^2 - 0.0137*temp + 10.62) - (2*sqrt(ionstrength))/
                    (1 + 1.4*sqrt(ionstrength)))

nonas$k1 <- with(nonas, 10^(-pk1))
nonas$k2 <- with(nonas, 10^(-pk2))

## calculate H+ concentration from measured pH
nonas$hplus <- with(nonas, 10^(-pH_ave))

## calculate a0-2 as per Kerri's gas flux calculations for 
##    proportions of DIC as CO2, HCO3 and CO3 respectively
nonas$ao <- with(nonas, 1/(1 + 10^(-pk1 + pH_ave) + 10^(-(pk1 + pk2) + 2*pH_ave)))
nonas$a1 <- with(nonas, 1/(10^(-pH_ave+pk1) + 1 + 10^(-pk2 + pH_ave)))
nonas$a2 <- with(nonas, 1/(10^(-2*pH_ave + (pk1 + pk2)) + 10^(-pH_ave + pk2) + 1))

## Null way of calculating total DIC: just add carb and bicarb, ignore CO2
## =================================================================================
## adding up carb and bicarb mg/L as Carbon
nonas$totaldicnull <- with(nonas, bicarbC + carbC)

## First way of calculating total DIC: use CO2 carbon calculated using k1
## =================================================================================
## following applied chemical equations from 
## http://www.soest.hawaii.edu/oceanography/courses/OCN623/Spring2012/CO2pH.pdf
##    assuming units are applicable
## k1 = [H+][HCO3-]/[CO2]
## k2 = [H+][CO3(2-)]/[HCO3-]
## The above have been rearranged to get the unknown (we know H+, k constants, and carb
##    and bicarb)
bicarbMOLG <- 1/(1.0079 + 12.0107 + 3*15.9994) # mol/g OR mmol/mg
nonas$co2mmol <- with(nonas, hplus*(Bicarb_mgL*bicarbMOLG)/k1) # mmol

nonas$co2mgL <- with(nonas, co2mmol*co2weight) # mg/L CO2
nonas$co2C <- with(nonas, co2mmol*12.0107) # mg/L C

with(nonas, co2mgL/(Bicarb_mgL + carb_mgL)) # proportion of co2 compared with the others

nonas$totaldic0 <- with(nonas, bicarbC + carbC + co2C)

## Second way of calculating total DIC: use a0-2 
## =================================================================================
# if a1 and a2 are proportions of HCO3 and CO3 of DIC then...
nonas$totaldic1 <- with(nonas, bicarbC/a1)
nonas$totaldic2 <- with(nonas, carbC/a2)

## Third way: combine a-values to use all data at once
## =================================================================================
nonas$totaldic3 <- with(nonas, (bicarbC+carbC)/(a1+a2))
## =================================================================================

## is there an outlier in carb_mgL ?
plot(density(nonas$carb_mgL, na.rm = TRUE))
# very high cond and Na and alk in this lake so might actually be legit

## plot relationships between total DIC estimates
## bonus, function that will plot a 1:1 line or other
panel.regression <- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                              cex = 1, col.regres = "red", ...) 
{ 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex) 
  ok <- is.finite(x) & is.finite(y) 
  if (any(ok)) 
    abline(0,1, col = col.regres, ...) # lm(y[ok] ~ x[ok]), col = col.regres, ...) 
} 
# scatterplot matrix of all options
pairs(nonas[,c('totaldicnull', 'totaldic0', 'totaldic1', 'totaldic2', 'totaldic3')]
      , panel=panel.regression)
# seems like totaldic2 is the most aberrant
# PS for our pH the null would be pretty much the same

## save data for CO2 fluxes
write.csv(nonas, '../data/private/phamdat-mod.csv', row.names = FALSE)

