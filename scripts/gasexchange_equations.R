## a collection of equations for gas exchange from EW's file
## chemical enhancement qu'appelle lakes 2008 to 2010.xlsx in /Data/Prairies

## given variables: temp ("temp"), pH ("ph"), pCO2[atm] ("kpa") 
## (370 as per Mauna Loa but could be checked if need to change recent values in time series), 
  ## bar prs (kP) ("bar"), alt (m) ("alt"), wind (m/s) ("wind"), (alk (mcEq/L) ("alk"), salinity ("salt").

## given constants:
r1a <- -2225.22 # all r1 are Khco3...
r1b <- -0.049
r1d <- 8.91
r1e <- 336.6

r2a <- 1346.24 # all r2 are kd
r2b <- -0.126
r2d <- -6.44
r2e <- -196.4

## equations

# =IF(G2<>" ",G2*7*0.0000025,AF$2*7*0.0000025)
ionstrength <- cond*7*0.0000025 

# =(0.000011*F2^2-0.012*F2+6.58)-(0.5*SQRT(K2))/(1+1.4*SQRT(K2))
pk1 <- (0.000011*temp^2 - 0.012*temp + 6.58) - (0.5*sqrt(ionstrength))/(1 + 1.4*sqrt(ionstrength))

# =(0.00009*F2^2-0.0137*F2+10.62)-(2*SQRT(K2))/(1+1.4*SQRT(K2))
pk2 <- (0.00009*temp^2 - 0.0137*temp + 10.62) - (2*sqrt(K2))/(1 + 1.4*sqrt(K2))

# =1/(1+10^(-L2+D2)+10^(-(L2+M2)+2*D2))
ao <- 1/(1 + 10^(-pk1 + ph) + 10^(-(pk1 + pk2) + 2*ph))

# =1/(10^(-D2+L2)+1+10^(-M2+D2))
a1 <- 1/(10^(-ph+pk1) + 1 + 10^(-pk2 + ph))

# =1/(10^(-2*D2+(L2+M2))+10^(-D2+M2)+1)
a2 <- 1/(10^(-2*ph + (pk1 + pk2)) + 10^(-ph + pk2) + 1)

# =1.11+0.016*F2-0.00007*F2^2
pkh <- 1.11 + 0.016*temp - 0.00007*temp^2

# =((O4+2*P4)*(R4)-1000000*10^(-D4)+1000000*10^(-14)+D4)
alk <- ((a1 + 2*a2)*(dic) - 1000000*10^(-ph) + 1000000*10^(-14) + ph)

# =1000000*(E2*0.000001-10^(-14+D2)+10^(-D2))/(O2+2*P2) ... mcM
# dic <- 1000000*(alk*0.000001 - 10^(-14 + ph) + 10^(-ph))/(a1 + 2*a2) in case at some point alk exists not dic. 
  # but this equation not tested to be correct yet. Further, Kerri did a different calculation for missing DIC 
  # values using conductivity-DIC relationship.... but I don't have this calculation at hand.... 
  # Kerri suggested using real data and original calculations whenever possible and then substituting whenever necessary

# =R2*N2 ... mcM; Kerri replaced this with CO2 vs pH relationship for places where dic etc missing or something, 
  # but advised me to use the original relationship!
co2 <- dic*ao

# =S2/(10^-Q2) .... mcAtm
pco2 <- co2/(10^-pkh)

# =IF(J2="",+(I2/101.325)*((+H2*0.000001)*10^(-Q2+6)),EXP(-0.12806*(+J2/1000))*((+H2*0.000001)*10^(-Q2+6))) ... equilibrium co2
co2eq <- exp(-0.12806*(+ alt/1000))*((+ kpa*0.000001)*10^(-pkh + 6))) 

# log(S2)
co2log <- log10(co2) 

# =(S2-U2)*0.8
flux <- (co2 - co2eq)*1.056 # !!! ok the value changes in the spreadsheet from top value of .8 to 1.056 
  # because kerri did calculations based on cole and caraco to determine the latter! 
  # reconsider if suddenly checking any other region's data e.g. boreal

# =Y2*X2
fluxenh <- flux*alpha

# =10^(-D2)
ah <- 10^(-ph)

# =2.07+0.215*F3^1.7 k cm/h .. no idea what this is
kspeed1 <- 2.07 + 0.215*wind^1.7

# =G2/60/60 k cm/s ..
kspeed2 <- kspeed1/60/60

# =10^-(0.000152*(C2^2)-0.0129*C2+6.5784)
k1 <- 10^-(0.000152*(temp^2) - 0.0129*temp + 6.5784)

# =10^-(0.000121*(C2^2)-0.0149*C2+10.626)
k2 <- 10^-(0.000121*(temp^2) - 0.0149*temp + 10.626)

# =10^-(0.0002*C2^2-0.0444*C2+14.953)
kw <- 10^-(0.0002*temp^2 - 0.0444*temp + 14.953)

# =(0.61563+0.05316*C2)*0.00001 D cm2/s .. no idea what D is...
dspeed <- (0.61563+0.05316*temp)*0.00001

# note that if salinity data exist, this will not need to run!!
salt <- cond*.68/1000

# =EXP($X$2+$Y$2*$T2^0.5+(10^4)*$Z$2/(C2+273)+$AA$2*LN(C2+273)) 
r1 <- exp(r1a + r1b*salt^0.5 + (10^4)*r1d/(temp + 273) + r1e*ln(temp + 273))

# =EXP($X$3+$Y$3*$T2^0.5+(10^4)*$Z$3/(C2+273)+$AA$3*LN(C2+273))
r2 <- exp(r2a + r2b*salt^0.5 + (10^4)*r2d/(temp + 273) + r2e*ln(temp + 273))

# =M2+N2*K2/E2
r <- r1 + r2*kw/ah

# =(O2*R2/L2)^0.5
q <- (r*t/dspeed)^0.5

# =L2/H2
z <- dspeed/kspeed2

# =1+E2^2*(I2*J2+I2*E2)^-1
t <- 1 + ah^2*(k1*k2 + k1*ah)^-1

# =R2/((R2-1)+TANH(P2*Q2)/(P2*Q2))
alpha <- t/((t - 1) + tanh(q*z)/(q*z))

## other comments: bars are plugged in in the calculations but it is a function of altitude.... 
  # and can't find it in the database!! So, can calculate in R or grab annual mean values from weather website as kerri suggests
# bar <- 1e^(-alt/7) # 7 is the "scale height" of the atmosphere, 1 is the pressure at sea level; alt needs to be in km!!

