## discovered that salinity values prior to 2003 seemed lower and more variable across the lakes
## --> started to try to figure out why
## YSI models through time: copy-paste from email to Peter: "I went in but could only 
##    find 2004 and 2007 out of the older ones [lab manuals] - which both used the YSI 85. 
##    Vince told me that in the late -90s (1998?) the probes had been changed, and that 
##    since 2013 we've had YSI pro plus."
## --> change from 2002 to 2003 was unlikely to be due to a change in hardware
## --> checked specs for YSI85 and discovered that salinity is not directly measured,
##    but calculated based on temperature and pressure
## --> emailed YSI to see if I could get the calculations that the meter does internally
## They emailed it to me and this file is saved in data/private right now. The calculations 
##    are copied into R in this script (checked for reproducibility)

## download necessary background data
params <- readRDS("data/private/params.rds")

## take out most likely outliers (see metadata file in shared google drive)
params <- subset(params, Salinity < 2 & pH < 12 & pH >= 7 & Conductivity < 2500)
removed <- which(params$Lake == "D" & params$Conductivity > 1000)
params <- params[-removed,]
removek <- which(params$Lake == "K" & params$Conductivity < 500)
params <- params[-removek,]
removel <- which(params$Lake == "L" & params$Conductivity < 750)
params <- params[-removel,]
removec <- which(params$Lake == "C" & params$Conductivity < 500)
params <- params[-removec,]


## load necessary packages
library("ggplot2")
library("reshape2")

## a few plots of measured data for the curious:
##    just salinity
  salplot <- ggplot(data = params, aes(x= Date, y = Salinity,
                                          group = Lake, colour = Lake)) +
    geom_line() +
    geom_point() +
    facet_wrap( "Lake" ) +
    theme(legend.position = "top")
  
  salplot

##    just conductivity
  condplot <- ggplot(data = params, aes(x= Date, y = Conductivity,
                                        group = Lake, colour = Lake)) +
    geom_line() +
    geom_point() +
    facet_wrap( "Lake" ) +
    theme(legend.position = "top")
  
  condplot
  
##    just pH
  phplot <- ggplot(data = params, aes(x= Date, y = pH,
                                      group = Lake, colour = Lake)) +
    geom_line() +
    geom_point() +
    facet_wrap( "Lake" ) +
    theme(legend.position = "top")
  
  phplot

  
##    just temperature
    tempplot <- ggplot(data = params, aes(x= Date, y = Temperature,
                                        group = Lake, colour = Lake)) +
    geom_line() +
    geom_point() +
    facet_wrap( "Lake" ) +
    theme(legend.position = "top")
  
  tempplot

##    conductivity and salinity plotted together
##    lookie here: https://github.com/hadley/ggplot2/wiki/Align-two-plots-on-a-page
  sal <- params[,c('Salinity', 'Date', 'Lake')]
  names(sal)[1:2] <- c("y", "Date")
  sal$panel <- rep("SAL")
  cond <- params[,c('Conductivity', 'Date', 'Lake')]
  names(cond)[1:2] <- c("y", "Date")
  cond$panel <- rep('COND')
  super <- rbind(sal, cond)
  
  p <- ggplot(data = super, mapping = aes(x = Date, y = y, colour = Lake)) # 
  p <- p + facet_grid(Lake+panel~., scale="free")
  p <- p + layer(data = sal, geom = c("point"), stat = "identity")
  p <- p + layer(data = cond, geom = c("point"), stat = "identity") + theme(legend.position = "none")
  p

## YSI85 calcs, which are based on conductivity and temperature readings on the probe:
## dbar: input pressure in decibars
## rratio: ratio of measured conductivity to the conductivity of 
##    the Standard Seawater Solution 
## temp: input deg C
## cond: input conductivity as mS/cm

## assuming 1m depth adds .1 bar to air pressure, can take rough mean of all sites 
##    and add that value for dbar; converted kpa to dbar
kpa <- mean(params$Pressure)
dbar <- kpa/10 + 1

##create function that spits out salinity
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

## call this function and add column to params
params <- transform(params, SalCalc = salcalc(Temperature, Conductivity, dbar))

## plot calculated salinity vs measured salinity
Index <- c("calculated", "measured")
salcalcplot <- ggplot(data = params, aes(x= Date, y = SalCalc, 
                                         group = Lake, colour = Index[1])) +
  ylab("Salinity (ppt)") +
  geom_line() +
  geom_point() +
  geom_line(data = params, aes(x = Date, y = Salinity, group = Lake, colour = Index[2])) +
  facet_wrap( "Lake" ) +
  theme(legend.position = "top")

salcalcplot

diffplot <- ggplot(data = params, aes(x= Date, y = SalCalc - Salinity, 
                                      group = Lake, colour = Lake)) +
  geom_abline() +
  geom_point( size = 1, shape= 21, fill="white") +
  facet_wrap( "Lake" ) + 
  ylab("Calculated - measured salinity (ppt)")

diffplot



## create pdfs to show:
pdf("data/private/salcalcplot.pdf", width = 15)
salcalcplot
dev.off()

pdf("data/private/salcalcdiffplot.pdf", width = 15)
diffplot
dev.off()


