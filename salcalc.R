## create function that spits out salinity, following conversations with YSI people and getting
##    the internal probe calculations for the variable --- these to replace all other salinity values!

salcalc <- function(temp, cond, dbar) {
  temp <- temp
  cond <- cond/1000 # in params, cond is in uS/cm, we need mS/cm
  dbar <- dbar
  consts <- read.csv("data/private/ysi85constants.csv", row.names = 1)
  consts <- as.data.frame(t(consts))
  attach(consts)
  
  # =B6/42.914
  rratio <- cond/42.914
  
  # =1+((B8*(N22+N23*B8+N24*B8^2))/(1+K22*B7+K23*B7^2+(K24+K25*B7)*B4))
  rp <- 1 + ((dbar*(e1 + e2*dbar + e3*dbar^2))/(1 + d1*temp + d2*temp^2 + 
                                                  (d3 + d4*temp)*rratio))
  # =H22+H23*B7+H24*B7^2+H25*B7^3+H26*B7^4
  rt1 <- c0 + c1*temp + c2*temp^2 + c3*temp^3 + c4*temp^4
  
  # =B4/(B10*B12)
  rt2 <- rratio/(rp*rt1)
  
  # =((B7-15)/(1+E28*(B7-15)))*(E22+E23*(B14^0.5)+E24*B14+E25*
  #   (B14^1.5)+E26*(B14^2)+E27*(B14^2.5))
  deltas <- ((temp - 15)/(1 + k*(temp - 15)))*(b0 + b1*(rt2^.5) + b2*rt2 + 
                                                 b3*(rt2^1.5) + b4*(rt2^2) +
                                                 b4*(rt2^2.5))
  
  # =B22+B23*(B14^0.5)+B24*B14+B25*(B14^1.5)+B26*(B14^2)+B27*(B14^2.5)+B17
  salinity <- a0 + a1*(rt2^.5) + a2*rt2 + a3*(rt2^1.5) + a4*(rt2^2) +
    a5*(rt2^2.5) + deltas
  
  detach(consts)
  salinity
}

