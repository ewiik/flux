## plot and analyse BP diel data

## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
    source("diel-buffalo.R")}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')

regvars <- readRDS('../data/private/regvars.rds')
params <- readRDS('../data/private/params-flux.rds')
buffreg <- subset(regvars, select = c('Date', 'Lake', 'Month', 'DOY', 'pH_surface',
                                      'lakepCO2','co2Flux'), Lake == 'B' & Year > 2013)

## load packages
library("ggplot2")

## cleaning dates in 2014 
dates <- c('12/08/2014', '15/07/2014', '19/07/2014', '26/06/2014') # last may not be a cleaning
#   date
dates <- as.POSIXct(dates, format = "%d/%m/%Y")
datesDOY <- format(dates, "%j")

## explore data by plotting
## plots
bdatf <- bdat
bdatf$Month <- as.factor(bdatf$Month)

ggplot(data=bdat, aes(y = co2corr, x = ODOrel1, col = ifelse(Hour >= 8 & Hour <=20, 
                                                                     'red', 'black'))) +
  scale_color_identity() +
  geom_point() #+
#geom_vline(xintercept = c(datesDOY))

ggplot(data=bdat, aes(y = co2corr, x = pco2, col = factor(Month))) 
# ifelse(Hour >= 8 & Hour <=20, 
#   'red', 'black')
#scale_color_identity() + # means it understands red and black in ifelse
#scale_color_manual(values=c('black', 'red')) +
geom_point() +
  geom_abline(intercept=0,slope=1) +
  ylab("pCO2 (measured)") +
  xlab("pCO2 (calculated)") +
  geom_text(label='2014', x=600, y=1600, col='black')

ggplot(data=bdat, aes(y = co2corr, x = DOY, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  geom_vline(xintercept = as.numeric(datesDOY)) +
  ylab("CO2 (ppm)") +
  xlab("Day of year 2014 (vertical lines cleaning dates)")

ggplot(data=bdatf, aes(y = ph1, x = Time, col = isDay)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  #geom_vline(xintercept = as.numeric(datesDOY)) +
  ylab("pH") +
  xlab("Hour")

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

## compare 2014 data with calculated data
real2014 <- subset(bdat, select = c('DOY', 'ODOrel1','co2corr', 'Hour', 'ph1'), 
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
recalcflux <- with(recalc, gasExchangeUser(temp = Temperature, cond = Conductivity, ph=pHalt, 
                                           wind=meanWindMS, alknotdic = FALSE,
                                           , dic = TICalt, kpa=Pressure,
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

## melt for various purposes
bdatmelt <- melt(bdat[,c("datetime", "ODOrel1", "ODOrel2")], id.vars = "datetime")
bdatmelt$variable <- as.character(bdatmelt$variable)
bdatmelt$variable[grep("1", bdatmelt$variable)] <- "1m"
bdatmelt$variable[grep("2", bdatmelt$variable)] <- "Deeper"
names(bdatmelt)[grep("var", names(bdatmelt))] <- "Depth"

oxyplot <- ggplot(data=bdatmelt, mapping = aes(y=value, x=datetime, group=Depth, lty=Depth)) +
  geom_line() +
  ylab("Oxygen saturation (%)") +
  xlab("Date")

## any DIC concentration evidence for carbonate precipitation?
ggplot(data=params[params$Lake =='B',], aes(y = TIC, x = Month, group = Year)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  facet_wrap('Year') +
  ylab("DIC (mg/L)") +
  xlab("Month")
