source("../functions/gasExchangeExtra.R")
## note on temp computer this is modified to take in alk !!
dictest <- 32.2/0.012
dictest2 <- 33.2/0.012
testing <- data.frame(temp = c(18,18), cond = c(800,800), ph = c(7.6, 7.9), dic=c(dictest2, dictest), 
                      alt=c(300, 300), salt=c(0,0), pco2atm=c(380, 380))
fluxes <- with(testing, gasExchangeExtra(kerri=FALSE, altnotkpa=TRUE, temp=temp, cond=cond, ph=ph,
                                        dic=dic,alt=alt, salt=salt,pco2atm = pco2atm, wind=4))

marls <- read.csv('../data/private/marls.csv')
marls$alt <- rep(0)
marls$alt[marls$lake == 'ct'] <- 138
marls$alt[marls$lake == 'ht'] <- 8
marls$alt[marls$lake == 'mt'] <- 375

marls$alkcarb <- marls$alk*0.6 # alk as mg/L caco3 in my data set --> co3
marls$alkmEq <- marls$alk/50
marls$alkuEq <- marls$alkmEq * 1000
## FIXME: not sure which unit of alk I need to calculate DIC using the formula in Flex
marlflux <- with(marls, gasExchangeExtra(temp = temp, cond = cond, ph=ph, wind=4,altnotkpa = TRUE, salt=0,
                                        alk = alkuEq, alt=alt, pco2atm = 380, alknotdic = TRUE))

marlsall <- cbind(marls, marlflux)
with(marlsall, plot(fluxenh ~ lake))
