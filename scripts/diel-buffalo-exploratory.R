## script just to plot and check BP diel data...
##    1. Looking at temporal trends
##    2. Comparing calculated vs measured CO2
##    3. Compaing rfu with measured routines chl (see also email from Sep 20 2016; Kimberly Gilmour)

## load packages
library("ggplot2")
library("gridExtra")
library("viridis")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  source("diel-buffalo.R")}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')

regvars <- readRDS('../data/private/regvars.rds')
params <- readRDS('../data/private/params-flux.rds')
buffreg <- subset(regvars, Lake == 'B' & Year > 2013)
buffchl <- subset(regvars, select = c('Date', 'Lake', 'Month', 'DOY', 'Chl_a_ug_L',
                                      'Chl_a_ug_L_sur'), Lake == 'B' & Year > 2013)
bdatf <- bdat
bdatf$Month <- as.factor(bdatf$Month)

## =========================================================================================
## 1. Looking at temporal trends, quality, cleaning, in the data
## =========================================================================================
## cleaning dates in 2014 
dates <- c('12/08/2014', '15/07/2014', '19/07/2014', '26/06/2014') 
# last may not be a cleaning date
dates <- as.POSIXct(dates, format = "%d/%m/%Y", tz="Canada/Saskatchewan")
datesDOY <- format(dates, "%j")

## a few temporal plots
ggplot(data=bdat, aes(y = co2corr, x = DOY, col = isDay)) + 
  scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  geom_vline(xintercept = as.numeric(datesDOY)) +
  ylab("CO2 (ppm)") +
  xlab("Day of year 2014 (vertical lines cleaning dates)")

ggplot(data=bdatf, aes(y = ph1, x = Time, col = isDay)) + 
  geom_point() +
  ylab("pH") +
  xlab("Hour")

## relationships between variables
ggplot(data=bdat, aes(y = co2corr, x = ODOrel1, col = ifelse(isDay, 'red', 'black'))) +
  scale_color_identity() +
  geom_point() #+

## differences between night and day
dielplot <- ggplot(data=bdat, aes(y = co2corr, x = ODOrel1, col=TimeofDay)) +
  geom_point() +
  scale_color_viridis(discrete=TRUE, alpha = 0.4, option='viridis', begin= 0.3, end=0.8, direction = -1) +
  #scale_color_distiller(palette = 'Greens') +
  ylab(expression(paste(italic(p)*"CO"[2]~"(ppm)"))) +
  xlab(expression(paste("O"[2]~"(%)"))) +
  geom_abline(slope=0,intercept = 400, linetype='dotted') +
  geom_vline(xintercept = 100, linetype='dotted') +
  theme_bw(base_size = 10) +
  theme(legend.title=element_blank(), legend.position='top')
ggsave(plot=dielplot, filename='../docs/private/bp-o2-co2.png')

## =============================================================================================
## 2. compare measured data with calculated data
## =============================================================================================
## all months and times of day
ggplot(data=bdatf, aes(y = co2corr, x = pco2, col = Month)) +
  geom_point(alpha=0.2) +
  geom_abline(intercept=0,slope=1) +
  ylab("pCO2 (measured)") +
  xlab("pCO2 (calculated)")

## subset diel data to match most likely times when lake sampled for routines,and 
##    calculate mean ph for this time interval
real2014 <- subset(bdat, select = c('DOY', 'ODOrel1','co2corr', 'Hour', 'ph1'), 
                   Hour >= 10 & Hour <= 13)
realsplit<- with(real2014[,c('DOY','ph1')], split(real2014[,c('DOY','ph1')], list(DOY)))
realmeans <- do.call(rbind, lapply(realsplit, colMeans, na.rm=TRUE))

## merge with routines 2014 data for same time interval
bp2014 <- subset(buffreg, select=c('DOY', 'lakepCO2','pH_surface'))
all2014 <- merge(bp2014, realmeans)

## subset routines data to same DOY
take <- which(bp2014$DOY %in% real2014$DOY)
bp2014sub <- bp2014[take,]

## is the difference just due to pH differences? replace routines pH with corresponding diel pH
paramssub <- params[params$Lake == 'B' & params$Year == 2014,]
take <- which(paramssub$DOY %in% all2014$DOY)
recalc <- paramssub[take,]

recalc$pHalt <- c(8.65, 8.96, 9.02, 8.73, 8.7)
recalc$TICalt <- 26.57 + 0.018 * recalc$Conductivity  # (mg/L)
recalc$TICalt <- recalc$TICalt / 0.012 # --> uM
recalcflux <- with(recalc, gasExchangeUser(temp = Temperature, cond = Conductivity, ph=pHalt, 
                                           wind=meanWindMS, alknotdic = FALSE,
                                           dic = TICalt, kpa=Pressure,
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
                      option = 'viridis', direction = -1) +
  ylab(expression(paste(mu*'atm'))) +
  xlab("Day of Year")
prec


## ===========================================================================================
## 3. does the chl extrapolated from rfu correspond to chl with our lab chl? + DOC/turb
## ===========================================================================================
chltest <- merge(buffchl, bdat, by='DOY') # we have five days
chltest2 <- subset(chltest, Hour >= 10 & Hour <=14)

chlsub <- data.frame(DOY=unique(chltest$DOY), Chl_a_ug_L = unique(chltest$Chl_a_ug_L), 
                     Chl_a_ug_L_sur = unique(chltest$Chl_a_ug_L_sur), x=11)
sur <- ggplot(chltest2, aes(y=Chl_a_ug_L_sur , x= chl, col=DOY)) +
  papertheme +
  geom_point() +
  scale_color_viridis() +
  ylab("Routines Surface Chl a") +
  theme(axis.title.x= element_blank(), legend.key.width=unit(2,"cm"))
  
int <- ggplot(chltest2, aes(y=Chl_a_ug_L, x= chl, col=DOY)) +
  papertheme +
  scale_color_viridis(guide=FALSE) +
  geom_point() +
  ylab("Routines Integrated Chl a") +
  xlab("Diel Chl from 10:00-14:00") 
grid.arrange(sur, int)

doctest <- merge(buffreg[,c('DOC_mg_L','DOY')], bdat, by='DOY') # we have five days
doctest2 <- subset(doctest, Hour >= 10 & Hour <=14)
ggplot(doctest2, aes(y=DOC_mg_L, x= cdom, col=DOY)) +
  papertheme +
  scale_color_viridis() +
  geom_point() +
  theme(legend.key.width=unit(2,"cm")) +
  ylab("Routines DOC (mg/L)") +
  xlab("Diel CDOM (ug/L) from 10:00-14:00") 
  
ggplot(doctest2, aes(y=DOC_mg_L, x= log10(turb), col=DOY)) +
  papertheme +
  scale_color_viridis() +
  geom_point() +
  theme(legend.key.width=unit(2,"cm")) +
  ylab("Routines DOC (mg/L)") +
  xlab("Diel turbidity (log10 NTU) from 10:00-14:00") 
