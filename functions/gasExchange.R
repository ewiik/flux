### functions for gas exchange in order

gasExchange <- function(temp, cond, ph, dic, alt, kpa, wind, salt) { 
    r1a <- -2225.22 # all r1 are Khco3...
    r1b <- -0.049
    r1d <- 8.91
    r1e <- 336.6

    r2a <- 1346.24 # all r2 are kd
    r2b <- -0.126
    r2d <- -6.44
    r2e <- -196.4
    localc <- 1.056
    
    ionstrength <- cond*7*0.0000025 

    pk1 <- (0.000011*temp^2 - 0.012*temp + 6.58) - (0.5*sqrt(ionstrength))/(1 + 1.4*sqrt(ionstrength))
    
    pk2 <- (0.00009*temp^2 - 0.0137*temp + 10.62) - (2*sqrt(ionstrength))/(1 + 1.4*sqrt(ionstrength))
    
    ao <- 1/(1 + 10^(-pk1 + ph) + 10^(-(pk1 + pk2) + 2*ph))
    
    a1 <- 1/(10^(-ph+pk1) + 1 + 10^(-pk2 + ph))
    
    a2 <- 1/(10^(-2*ph + (pk1 + pk2)) + 10^(-ph + pk2) + 1)
    
    pkh <- 1.11 + 0.016*temp - 0.00007*temp^2
    
    alk <- ((a1 + 2*a2)*(dic) - 1000000*10^(-ph) + 1000000*10^(-14) + ph) # all good up to here
    
    co2 <- dic*ao
    # kerri did a calculation for co2 based on pH which we need for patching NAs
    
    pco2 <- co2/(10^-pkh) # this also differs but is correct and can use that, kerri 
    # used a different relationship
    
    co2eq <- exp(-0.12806*(+ alt/1000))*((+ kpa*0.000001)*10^(-pkh + 6)) 
    
    co2log <- log10(co2) 
    
    ah <- 10^(-ph)
    
    kspeed1 <- 2.07 + 0.215*wind^1.7
    
    kspeed2 <- kspeed1/60/60
    
    k1 <- 10^-(0.000152*(temp^2) - 0.0129*temp + 6.5784)
    
    k2 <- 10^-(0.000121*(temp^2) - 0.0149*temp + 10.626)
    
    kw <- 10^-(0.0002*temp^2 - 0.0444*temp + 14.953)
    
    dspeed <- (0.61563+0.05316*temp)*0.00001
    
    if (missing(salt)) {
      salt <- cond*.68/1000 
      print("calculating salt")
    }
    
    r1 <- exp(r1a + r1b*salt^0.5 + (10^4)*r1d/(temp + 273) + r1e*log(temp + 273)) 
    
    r2 <- exp(r2a + r2b*salt^0.5 + (10^4)*r2d/(temp + 273) + r2e*log(temp + 273)) 
    
    r <- r1 + r2*kw/ah 
    
    t <- 1 + ah^2*(k1*k2 + k1*ah)^-1
    
    q <- (r*t/dspeed)^0.5
    
    z <- dspeed/kspeed2
    
    alpha <- t/((t - 1) + tanh(q*z)/(q*z))
    
    flux <- (co2 - co2eq)*localc # requires calculation for local conditions; one stored by EW is
    # good for prairies, sk
    
    fluxenh <- flux*alpha
    fluxenh
}

