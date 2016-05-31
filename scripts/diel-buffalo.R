## load in Helen's 2015 buoy data -- from file sent 28.01.2016
## And 2014 data from email sent Monday - April 25, 2016 
## two sonde depths, n.1 is the shallower one based on PAR..
## Not sure what cdom and mv are measures of.. and why so many temps?
## FIXME: really need to learn what's what with the column name matching (2015)
## FIXME: removed outliers but also 0s that are not real 0s.. need to remove (2015)

## load necessary packages
library("ggplot2")
library('StreamMetabolism')
library('mgcv')
library('viridis')
library('reshape2')

## load necessary functions
source('../functions/gasExchangeUser.R')
source('../data/private/salcalc.R')

## load in supplementary data
regvars <- readRDS('../data/private/regvars.rds')
params <- readRDS('../data/private/params-flux.rds')
buffreg <- subset(regvars, select = c('Date', 'Lake', 'Month', 'DOY', 'pH_surface',
                                      'lakepCO2','co2Flux'), Lake == 'B' & Year > 2013)
## mauna loa for 2014 months, such few data points I copy pasted manually from
##    ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt
ml <- read.csv('../data/maunaloa2014summer.csv')
names(ml)[which(names(ml)=='pCO2')] <- 'maunaloa'
ml <- ml[,-which(names(ml) == 'pCO2interp')]

## ======================================================================================
## 2014
## ======================================================================================

## read in data; initially got only part of the supporting info, but the co2 that Helen
##    already corrected. Then asked for pH etc and she sent that but it doesn't have
##    the corrected values
bdat2014 <- read.csv('../data/private/BuffPd_DataBlobs_Helen.csv', skip = 4, 
                     col.names = c('datetime', "winddir", "windsp", "airtemp", 
                                   "relhum", "pressure", "dailyrain", "par1", "par2",
                                   "par3", "temp1", "ODOrel1", "ODOabs1", "temp2",
                                   "temp3", "temp4", "temp5", "co2corr"))
bdat2014supp <- read.csv("../data/private/BPBuoy2014.csv", skip = 4, 
                         col.names = c('datetime',"batvolt", "winddir", "windsp", "airtemp", 
                                       "relhum", "pressure", "dailyrain", "par1", "par2",
                                       "par3", "cdom", "mvolts", "co2-1", "co2-2", "temp1",
                                       "cond1", "ph1", "ph1mV", "turb", "chl", "chlrfu", "bga1cell",
                                       "bga1rfu", "ODOrel1", "ODOabs1", "temp2", "cond2", "ph2", 
                                       "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3", 
                                       "temp4", "temp5", "temp6"))

## make time understandable
bdat2014 <- transform(bdat2014, datetime = as.POSIXct(as.character(datetime), 
                                              format = "%m/%d/%Y %H:%M"))
bdat2014 <- transform(bdat2014, Hour = as.numeric(format(datetime, format = "%H")))
bdat2014 <- transform(bdat2014, Month = as.numeric(format(datetime, format = "%m")))
bdat2014 <- transform(bdat2014, Day = as.numeric(format(datetime, format = "%d")))
bdat2014 <- transform(bdat2014, DOY = as.numeric(format(datetime, format = "%j")))
bdat2014 <- transform(bdat2014, Time = strftime(datetime, format = "%H:%M", tz="GMT-1"))
bdat2014 <- transform(bdat2014, isDay = ifelse(Hour < 21 & Hour > 7, TRUE, FALSE))

bdat2014supp <- transform(bdat2014supp, datetime = as.POSIXct(as.character(datetime), 
                                                      format = "%m/%d/%Y %H:%M"))
bdat2014supp <- transform(bdat2014supp, Hour = as.numeric(format(datetime, format = "%H")))
bdat2014supp <- transform(bdat2014supp, Month = as.numeric(format(datetime, format = "%m")))
bdat2014supp <- transform(bdat2014supp, Day = as.numeric(format(datetime, format = "%d")))
bdat2014supp <- transform(bdat2014supp, DOY = as.numeric(format(datetime, format = "%j")))
bdat2014supp <- transform(bdat2014supp, Time = as.numeric(format(datetime, "%H")) +
  as.numeric(format(datetime, "%M"))/60)
bdat2014supp <- transform(bdat2014supp, Week = as.numeric(format(datetime, "%U")))
bdat2014supp <- transform(bdat2014supp, isDay = ifelse(Hour < 21 & Hour > 7, TRUE, FALSE))

## in 2014 pressure is in HPA but need KPA for CO2 calcs
bdat2014supp$pressureKPA <- bdat2014supp$pressure/10

## cleaning dates in 2014 (sensor not cleaned in 2015)
dates <- c('12/08/2014', '15/07/2014', '19/07/2014', '26/06/2014') # last may not be a cleaning
#   date
dates <- as.POSIXct(dates, format = "%d/%m/%Y")
datesDOY <- format(dates, "%j")

## homogenise data
bdat2014supp$pressureKPAwater <- bdat2014supp$pressureKPA + .08
bdat2014supp <- merge(bdat2014supp, ml)

## calculate salinity
bdat2014supp <- transform(bdat2014supp, salcalc = salcalc(temp=temp1, cond=cond1, 
                                                         dbar=pressureKPAwater/10))
bdat2014suppflux <- with(bdat2014supp, gasExchangeUser(temp = temp1, cond = cond1, ph=ph1, 
                                                       wind=windsp, alknotdic = FALSE, 
                                                       salt = salcalc, kpa=pressureKPA,
                                                       pco2atm = maunaloa, diccalc = TRUE))
bdat2014suppall <- cbind(bdat2014supp, bdat2014suppflux)
bdat2014corr <- subset(bdat2014, select=c('datetime', 'co2corr'))
bdat2014full <- merge(bdat2014suppall, bdat2014corr, all.x = TRUE)

## let's look at when sun up and down
sundat <- sunrise.set(lat = 50.648016, long=-105.5072930, date = '2014/06/01', 
                      timezone = 'UTC+6', num.days=90)
sundat <- transform(sundat, Month = as.numeric(format(sunrise, format = "%m")))
sundat <- transform(sundat, Day = as.numeric(format(sunrise, format = "%d")))
sundat <- transform(sundat, UpTime = as.numeric(sundat$sunrise - trunc(sundat$sunrise, "days")))
sundat <- transform(sundat, DownTime = as.numeric(sundat$sunset - trunc(sundat$sunset, "days")))

## merge sundat with BP data
bdat2014full <- merge(bdat2014full, sundat)

## change isday to depend on sunup and down rather than rigid times
bdat2014full$isDay <- ifelse(bdat2014full$Time < bdat2014full$DownTime & 
                               bdat2014full$Time > bdat2014full$UpTime, 
                             TRUE, FALSE)

## save object for other purposes
saveRDS(bdat2014full, '../data/private/bpbuoy2014-mod.rds')

## plots
bdat2014fullf <- bdat2014full
bdat2014fullf$Month <- as.factor(bdat2014fullf$Month)

ggplot(data=bdat2014full, aes(y = co2corr, x = ODOrel1, col = ifelse(Hour >= 8 & Hour <=20, 
                                                                 'red', 'black'))) +
  scale_color_identity() +
  geom_point() #+
  #geom_vline(xintercept = c(datesDOY))

ggplot(data=bdat2014full, aes(y = co2corr, x = pco2, col = factor(Month))) 
  # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  geom_abline(intercept=0,slope=1) +
  ylab("pCO2 (measured)") +
  xlab("pCO2 (calculated)") +
  geom_text(label='2014', x=600, y=1600, col='black')

ggplot(data=bdat2014full, aes(y = co2corr, x = DOY, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
    #   'red', 'black')
    #scale_color_identity() + # means it understands red and black in ifelse
    scale_color_manual(values=c('black', 'red')) +
    geom_point() +
    geom_vline(xintercept = as.numeric(datesDOY)) +
    ylab("CO2 (ppm)") +
    xlab("Day of year 2014 (vertical lines cleaning dates)")

ggplot(data=bdat2014fullf, aes(y = ph1, x = Time, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  #geom_vline(xintercept = as.numeric(datesDOY)) +
  ylab("pH") +
  xlab("Hour")

ggplot(data=bdat2014full, aes(y = co2corr, x = ODOrel1, col=isDay)) +
  geom_point() +
  ylab('CO2 (ppm)') +
  xlab('Oxygen (%)') +
  geom_text(label='2014', x=8.8, y=1600, col='black')

## compare 2014 data with calculated data
real2014 <- subset(bdat2014full, select = c('DOY', 'ODOrel1','co2corr', 'Hour', 'ph1'), 
                   Hour >= 10 & Hour <= 15)
bp2014 <- subset(buffreg, select=c('DOY', 'lakepCO2','pH_surface'))
all2014 <- merge(bp2014, real2014)

take <- which(bp2014$DOY %in% real2014$DOY)
bp2014sub <- bp2014[take,]

with(all2014, plot(co2corr ~ DOY, ylab = 'pCO2 (ppm|uatm)', xlab = 'Day of year'))
points(all2014$lakepCO2 ~ all2014$DOY, pch=4)
abline(v = c(datesDOY))
legend('topright', legend = 'o=measured \n x=calculated', border = FALSE)

with(all2014, plot(ph1 ~ DOY, ylab = 'pH probe/routines', xlab = 'Day of Year', ylim = c(7.5,10)))
points(all2014$pH_surface ~ all2014$DOY, pch=4)
abline(v = c(datesDOY))
legend('topright', legend = 'o=measured \n x=calculated', border = FALSE)

## is the difference just due to pH differences?
paramssub <- params[params$Lake == 'B' & params$Year == 2014,]
take <- which(paramssub$DOY %in% all2014$DOY)
recalc <- paramssub[take,]

recalc$pHalt <- c(8.65, 8.96, 9.02, 8.73, 8.7)
recalc$TICalt <- 26.57 + 0.018 * recalc$Conductivity  # (mg/L)
recalc$TICalt <- recalc$TICalt / 0.012 # --> uM
recalcflux <- with(recalc, gasExchangeExtra(temp = Temperature, cond = Conductivity, ph=pHalt, 
                                           wind=meanWindMS, alknotdic = FALSE,
                                          kerri = FALSE, dic = TICalt, kpa=Pressure,
                                          pco2atm=397, salt=SalCalc))
recalc <- cbind(recalc, recalcflux)

## let's put all the pco2 for 2014 together for a ggplot
names(recalc)[which(names(recalc) == 'pco2')] <- 'Recalculated'
names(all2014)[which(names(all2014) %in% c('lakepCO2', 'co2corr'))] <-
  c('Calculated', 'Measured')
names(bp2014sub)[which(names(bp2014sub)== 'lakepCO2')] <-
  'Calculated'

bp2014m <- melt(bp2014sub, id.vars = 'DOY', measure.vars = 'Calculated',
                variable.name = 'whichCO2',
                factorsAsStrings = FALSE)
all2014m <- melt(all2014, id.vars = c('DOY'), 
                 measure.vars = c('Measured'),
                 variable.name = 'whichCO2',
                 factorsAsStrings = FALSE)
recalcm <- melt(recalc, id.vars = c('DOY'), 
                 measure.vars = 'Recalculated',
                 variable.name = 'whichCO2',
                 factorsAsStrings = FALSE)
recalcs <- rbind(recalcm, all2014m, bp2014m)

recalcs$whichCO2 <- factor(recalcs$whichCO2, levels = c("Measured", "Calculated","Recalculated"))

prec <- ggplot(data=recalcs, aes(y = value, x = DOY, colour = whichCO2)) + 
  theme_bw() +
  geom_jitter(width = 5) +
  scale_color_viridis(expression(paste(italic('p')*'CO'[2])), discrete = TRUE, 
                      option = 'plasma', direction = -1) +
  ylab(expression(paste(mu*'atm'))) +
  xlab("Day of Year")

## difference between max value and recalc value?
spldf <- with(all2014, split(all2014, list(DOY)))
spldf <- lapply(spldf, '[[', 'co2corr')
lapply(spldf, max)
recalcflux$pco2

## stats
lm2014night <- with(bdat2014full[bdat2014full$isDay == FALSE,], lm(co2corr~pco2))
lm2014day <- with(bdat2014full[bdat2014full$isDay == TRUE,], lm(co2corr~pco2))
lm2014int <- with(bdat2014full, lm(co2corr~pco2))
with(bdat2014full, cor(x=co2corr, y=pco2, use='complete.obs', method='pearson'))

# just with the 10am-3pm data
ten3 <- subset(bdat2014full, Hour >= 10 & Hour <=15)

## could do this to make converge but also did manually, required four iterations
#reps <- 8 #takes a while for it to converge
#for (i in 1:reps) {
#   start <- if (i == 1 ) {coef(glm(co2corr ~ pco2, data = ten3, family = gaussian))}
#     else { coef(glmmod) }
#   glmmod <- glm(co2corr ~ pco2, data = ten3, family = Gamma(link = "identity"), start = start,
#                 control = glm.control(maxit=100))
# }

start <- c(-69.0709700, 0.9070629 )
glmmod <- glm(co2corr ~ pco2, data = ten3, family = Gamma(link = "identity"), start = start,
                               control = glm.control(maxit=100))

pglm <- ggplot(data=ten3, aes(y = co2corr, x = pco2)) + 
  theme_bw() +
  geom_point() +
  geom_abline(slope=1, intercept=0, color='red', lty=3) +
  geom_ribbon() + ## FIXME: add the regression line here
  geom_abline(intercept = glmmod$coefficients[1], slope=glmmod$coefficients[2], col = 'black') +
  scale_color_viridis(discrete = TRUE, 
                      option = 'plasma', direction = -1) +
  ylab(expression(paste('Measured'~italic('p')*'CO'[2]~'('*mu*'atm)'))) +
  xlab(expression(paste('Calculated'~italic('p')*'CO'[2]~'('*mu*'atm)')))


# residual plot etc
plot(resid(lm2014int) ~ lm2014int$fitted.values, ylab='Residuals (pCO2 uatm)', 
     xlab='Fitted values (pCO2 uatm)')

with(bdat2014full, plot(co2corr ~ pco2, col=ifelse(isDay, 'red', 'black'),
                        xlab = 'Calculated pCO2', ylab = 'Measured pCO2'))
abline(0,1, lty=3)
abline(lm2014day$coefficients[1], lm2014day$coefficients[2], col = 'red')
abline(lm2014night$coefficients[1], lm2014night$coefficients[2], col = 'black')
#abline(-71.6, 1.143408, lty = 3, col = 'green')
legend('topleft', col = c('red', 'black', 'black'), pch=c(1,1, NA), lty = c(1,1,3), 
legend = c('daytime values and regression line', 
           'night-time values and regression line',
           '1:1 line'))
legend('bottomright', col = c('red', 'black'), title = 'R2', 
       legend = c('day = 0.72', 'night = 0.92'))

with(bdat2014full, plot(pco2 - co2corr ~ co2corr, xlab = 'Measured pCO2', 
                        ylab = 'Calculated - measured pCO2'))
abline(0,0)

## any DIC concentration evidence for carbonate precipitation?
ggplot(data=params[params$Lake =='B',], aes(y = TIC, x = Month, group = Year)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  facet_wrap('Year') +
  ylab("DIC (mg/L)") +
  xlab("Month")

## gamming the cyclic stuff in the data

m1 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(Week), data = bdat2014full)

m2 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(Week) + ti(Time, Week, bs = c("cc","tp")), 
          data = bdat2014full)
anova(m1, m2, test = "LRT")

m3 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")), 
          data = bdat2014full)
mco2 <- gam(co2corr ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")), 
          data = bdat2014full)

plot(m2, scheme=2, ylab='mean-0 pH', main='mean-0 pH') # latter works for 3D plot
#   former for the 2D plot
plot(acf(resid(m2)))

pdf("diel-bpgam.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(m2, pages = 4, scheme = 2)
par(op)
dev.off()

## predict for periods of varying interactions.
N <- 200
simDOY <- c(171:173, 211:213, 237:239)
DOYgroup <- factor(rep(1:3, times=c(3,3,3)))
DOYgroup <- factor(rep(c('171:173', '211:213', '237:239'), times=c(3,3,3)))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(bdat2014full$Time, na.rm=TRUE),
                                      max(bdat2014full$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
m3pred <- predict(m3, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, m3pred, DOYgroups)
names(predicted)[which(names(predicted)=='m3pred')] <- 'pH'
ggplot(predicted, aes(x = Time, y = pH, group= DOY, col=DOYgroups)) +
  theme_bw() +
  geom_line() +
  scale_colour_discrete(name="Day of Year") +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab('pH')

## for co2
N <- 200
simDOY <- c(168:170, 190:192, 228:230)
DOYgroup <- factor(rep(1:3, times=c(3,3,3)))
DOYgroup <- factor(rep(c('Jun 17-19', 'Jul 9-11', 'Aug 16-18'), times=c(3,3,3)))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(bdat2014full$Time, na.rm=TRUE),
                                      max(bdat2014full$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
mco2pred <- predict(mco2, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, mco2pred, DOYgroups)
names(predicted)[which(names(predicted)=='mco2pred')] <- 'pCO2'

labdatco2 <- data.frame(x = 3, y = 380, label = "Atmospheric mean")

co2modplot <- ggplot(predicted, aes(x = Time, y = pCO2, group= DOY, col=DOYgroups)) +
  theme_bw() +
  scale_colour_brewer(name = "Date", type = 'qual', palette = 'Dark2', direction=1) +
  geom_line() +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab(expression(paste(italic(p)*"CO"[2]~"(ppm)"))) +
  geom_abline(intercept = 398, slope = 0, linetype='dotted') +
  geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
            show.legend = FALSE, inherit.aes = FALSE)
ggsave('../docs/private/bp-diel-co2gam.png', width=20, height=15, units='cm')

## gam on other params
bdat2014fullf <- transform(bdat2014fullf, daylength = sunset-sunrise)
bdat2014fullf$daylength <- as.numeric(bdat2014fullf$daylength)
bpmod <- gam(ph1 ~
                    s(chl) +
                    s(airtemp) +
                    s(daylength) +
                    s(ODOrel1),
                  data = bdat2014fullf,
                  select = TRUE, method = "REML", family = gaussian(),
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, trace = TRUE,
                                        newton = list(maxHalf = 60)))
saveRDS(bpmod, '../data/private/bpmod.rds')

bpco2phmod <- gam(co2corr ~ s(ph1), method = "REML", family = gaussian,
                  na.action = na.exclude, data=bdat2014fullf,
                  control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
saveRDS(bpco2phmod, '../data/private/bpco2phmod.rds')
## ==========================================================================================
## 2015
## ==========================================================================================
bdat <- read.csv("../data/private/BPBuoyData2015raw.csv", skip = 4, 
                 col.names = c("datetime", "batvolt", "winddir", "windsp", "airtemp", 
                               "relhum", "pressure", "dailyrain", "par1", "par2",
                               "par3", "cdom", "mvolts", "co2-1", "co2-2", "temp1",
                               "cond1", "ph1", "ph1mV", "turb", "chl", "chlrfu", "bga1cell",
                               "bga1rfu", "ODOrel1", "ODOabs1", "temp2", "cond2", "ph2", 
                               "ph2mV", "bga2cell", "bga2rfu", "ODOrel2", "ODOabs2", "temp3", 
                               "temp4", "temp5", "temp6", "temp7"))

press <- read.csv("../data/bp-pressure2015.csv")
names(press)[which(names(press) == 'Pressure')] <- 'AirPressure'
press <- press[,-which(names(press) == 'Year')]

bdat <- transform(bdat, datetime = as.POSIXct(as.character(datetime), 
                                              format = "%m/%d/%Y %I:%M:%S %p"))
bdat <- transform(bdat, Hour = as.numeric(format(datetime, format = "%H")))
bdat <- transform(bdat, Month = as.numeric(format(datetime, format = "%m")))
bdat <- transform(bdat, Day = as.numeric(format(datetime, format = "%d")))
bdat <- transform(bdat, DOY = as.numeric(format(datetime, format = "%j")))

bdat <- transform(bdat, isDay = ifelse(Hour < 21 & Hour > 7, TRUE, FALSE))

## get pressure in
bdat <- merge(bdat, press, by=c('Month', 'Day'))

## check the CO2 columns
with(bdat, plot(co2.1 ~ datetime))

## not sure why all june and july values are completely different for both co2 profiles..!?     
with(bdat, plot(density(co2.1), main='Density of 2015 uncorrected all-outliers-in pCO2'))
with(bdat, plot(density(co2.2)))

## remove the crazy negative and positive values and bound the ppm
##    but first make note of which rows these are
outliers <- which(abs(bdat$co2.1) > 1200 | abs(bdat$co2.2) > 1200)
bdat<- bdat[-outliers,]
## still bimodal but don't wanna remove too much in case it's real

## initial pressure in hg inches or kPa (met), need to get hPa for corrections
bdat$pressure[grep('-Invalid-', bdat$pressure)] <- NA
bdat$pressure <- as.character(bdat$pressure)
bdat$pressure <- as.numeric(bdat$pressure)
bdat$pressureHPA <- bdat$pressure/0.0295299830714

bdat$AirPressureHPA <- bdat$AirPressure * 10

## deal with similar issues with wind (is in m/s)
bdat$windsp[grep('-Invalid-', bdat$windsp)] <- NA
bdat$windsp <- as.character(bdat$windsp)
bdat$windsp <- as.numeric(bdat$windsp)

## plot some diagnostics and look if pH and oxygen determining CO2
##    assuming I know which odo and ph measure relates to which co2 probe...
with(bdat, plot(co2.1 ~ Hour, col = ifelse(isDay, "red", "black")))
with(bdat, plot(co2.2 ~ Hour, col = ifelse(isDay, "red", "black")))

with(bdat, plot(co2.1 ~ ph1, col = ifelse(Month < 9, "red", "black")))
with(bdat, plot(co2.1 ~ ODOrel1, col = ifelse(isDay, "red", "black")))
## something wrong with ODO measurements!! loads 0s. Let's remove 0s
## Existing data mostly contiguous in time
with(bdat[-which(bdat$ODOrel1 == 0), ], plot(co2.1 ~ ODOrel1, col = ifelse(Day, "red", "black")))


with(bdat, plot(co2.2 ~ ph2, col = ifelse(Month < 9, "red", "black")))
with(bdat, plot(co2.2 ~ ODOrel2, col = ifelse(Month < 9, "red", "black")))
## here, something fishy with the month of September.. CO2 mostly 0

with(bdat[-which(bdat$ODOrel1 == 0), ], 
     plot(ODOrel1 ~ ph1, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
with(bdat, plot(ODOrel2 ~ ph2, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
## night and day different slopes in sonde1, doesn't happen in sonde2 even when
##    following the ODOrel1 subset

with(bdat[-which(bdat$ODOrel1 == 0), ], 
     plot(ODOrel1 ~ chl, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))
with(bdat, plot(ph1 ~ chl, col = ifelse(Hour > 21 | Hour < 11, "black", "red")))


## correction as per ecohydrology paper
#The post-measurement correction of sensor output as pCO2 related to changes in temperature 
#and pressure were determined empirically for the Vaisala GMT IRGAs (Vaisala Oyj, 2008). 
#Pre-corrected sensor output needs to be reduced by 0.15% of the measured reading per hPa 
#increase in pressure relative to calibration pressure (typically 1013 hPa). Pressure readings
#below the calibration pressure require increasing the sensor output by 0.15% of the measured
#reading per hPa. Pre-corrected sensor output also needs to be increased by 0.3% of the measured 
#reading per ?C of increased temperature relative to calibration temperature (typically 25?C). 
#Water temperature above the calibration temperature requires decreasing the sensor output by 
#0.3% of the measured reading per ?C. An additional correction is required when sensors are
#deployed in an aquatic environment. Water depth affects the pressure exerted on the sensor, 
#and this water depth correction is added to the atmospheric pressure correction. For example, 
#10 cm water depth above the sensor location corresponds to an increase in pressure of 9.81 hPa. 
#If sensor output is 1000 ppm pCO2, the output needs to be reduced by 14.72 ppm because of 
#the additional pressure exerted by the increased water level. If the water depth varies 
#relative to the sensor location, a depth correction must be determined for each measurement 
#time interval (e.g. by recording water depth using a pressure transducer or other water level sensor)
corrpress <- bdat$co2.1*((1013 - bdat$AirPressureHPA) * 0.0015)
corrtemp <- bdat$co2.1*((25 - bdat$temp1) * 0.003)
corrairpress <- -bdat$co2.1*((1014.14 - 1013.25) * 0.0015)
# 89cm pressure calculated at http://www.calctool.org/CALC/other/games/depth_press
#   is where the CO2 sensor is.

bdat$co2.1corr <- bdat$co2.1 + corrpress + corrtemp + corrairpress

## plot pH with co2 
bdat2 <- bdat
bdat2$Month <- as.factor(bdat2$Month)
ggplot(data=bdat2, aes(y = co2.1corr, x = ph1, col = Month)) +
  geom_point() +
  ylab('CO2 (ppm) corrected') +
  xlab('pH') +
  geom_text(label='2015', x=8.9, y=1250, col='black')

## calculate CO2 for 2015
bdat$dic <- rep(NA)
bdatflux <- with(bdat, gasExchangeExtra(temp = temp1, cond = cond1, ph=ph1, wind=windsp, alknotdic = FALSE,
                                       kerri = TRUE, dic = dic, kpa=AirPressure))
bdatall <- cbind(bdat, bdatflux)

ggplot(data=bdatall, aes(y = co2.1corr, x = pco2, col = factor(Month))) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  geom_abline(intercept=0,slope=1) +
  ylab("pCO2 (measured)") +
  xlab("pCO2 (calculated)") +
  geom_text(label='2015', x=200, y=1250, col='black')

lm2015 <- with(bdatall, lm(co2.1corr~pco2 -1))
with(bdatall, cor(co2.1corr, pco2,use='complete.obs', method='pearson'))

## plot 2015 data
ggplot(data=bdat, aes(y = co2.1, x = DOY, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  ylab("CO2 (ppm)") +
  xlab("Day of year 2015")

