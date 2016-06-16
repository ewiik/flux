### Get site specific data/predictions

## load packages
library('ggplot2')
library('viridis')
library('cowplot')
library('mgcv')
library('extrafont')
library('gridExtra')

## Set defaults
theme_set(theme_bw())

## load in data
regvars <- readRDS("../data/private/regvars.rds")

## Load in gam models
co2mod <- readRDS("../data/private/co2mod.rds")
egmod.red2 <- readRDS("../data/private/egmodred2.rds")
co2modnull <- readRDS("../data/private/co2modnull.rds")

## change to match what was used to create the models
regvarf <- regvars
regvarf <- transform(regvarf, Year = as.factor(Year)) # make Year into factor for re
regvarf2 <- regvarf
regvarf2$`Chl_a_ug_L`[regvarf2$`Chl_a_ug_L` <=0 ] <- 0
regvarf2$`DOC_mg_L`[regvarf2$`DOC_mg_L` <=0 ] <- 0
regvarf2$`Chl_a_ug_L` <- regvarf2$`Chl_a_ug_L` + 1
regvarf2$`DOC_mg_L` <- regvarf2$`DOC_mg_L` + 1
regvarf2 <- transform(regvarf2, dummy = rep(1, nrow(regvarf2)))

## what is mean pH? perhaps add a dotted line at this point in the plots
meanpH <- with(regvarf2[which(regvarf2$pH_surface >7 & regvarf$pH_surface < 11.1),], 
               mean(pH_surface, na.rm=TRUE)) # some outliers in there
meanco2 <- with(regvarf2, 
               mean(co2Flux, na.rm=TRUE))
meanO2 <- with(regvarf2, 
                mean(Oxygen_ppm, na.rm=TRUE))
meanchl <- with(regvarf2, 
                mean(Chl_a_ug_L, na.rm=TRUE))
meanTDN <- with(regvarf2, 
                mean(TDN_ug_L, na.rm=TRUE))
meanGPP <- with(regvarf2, 
                mean(GPP_h, na.rm=TRUE))

## create a theme to save linespace in plots
papertheme <- theme_bw(base_family='Arial') +
  theme(legend.position='top') 

## generate predicted for co2 and ph for all signif variables
## co2 ~ ph: full model
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

co2.pdat <- with(droplevels(regvarf2),
                  data.frame(`pH_surface` = rep(seq(min(`pH_surface`, na.rm = TRUE),
                                                    max(`pH_surface`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
co2.pdat <- merge(co2.pdat, lakeXbar)
co2.pred <- predict(co2mod, newdata = co2.pdat, type = "terms", se.fit = TRUE)

whichCols <- grep("pH", colnames(co2.pred$fit))
whichColsSE <- grep("pH", colnames(co2.pred$se.fit))
co2.pdat <- cbind(co2.pdat, Fitted = rowSums(co2.pred$fit[, whichCols]), 
                  se.Fitted = rowSums(co2.pred$se.fit[, whichColsSE]))
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

co2.pdat <- with(co2.pdat, transform(co2.pdat, Fittedplus = Fitted + se.Fitted))
co2.pdat <- with(co2.pdat, transform(co2.pdat, Fittedminus = Fitted - se.Fitted))

## make into original limits
shiftco2 <- attr(predict(co2mod, newdata = co2.pdat, type = "iterms"), "constant")
co2.pdatnorm <- co2.pdat
co2.pdatnorm <- with(co2.pdatnorm, transform(co2.pdatnorm, Fitted = Fitted + shiftco2))
labdatco2 <- data.frame(x = 7.5, y = meanco2 - 18, label = "mean flux")

## add quantiles
phquants <- quantile(regvarf2$pH_surface, c(.05,.95), na.rm = TRUE)

co2plot <- ggplot(co2.pdatnorm, aes(x = pH_surface, y = Fitted, 
                    colour = ifelse(Lake == "L", "Last Mountain", 
                                    ifelse(Lake == "B","Buffalo Pound",
                                           ifelse(Lake == "C", "Crooked",
                                                  ifelse(Lake == "D", 
                                      "Diefenbaker", "Wascana, Pasqua,\nKatepwa")))))) +
  annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = .2) +
  geom_line() +
  theme_bw(base_size = 12, base_family = 'Arial') +
  #geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
   #           alpha = 0.25, fill='white') +  
  scale_colour_brewer(name = "Lake", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE, inherit.aes = FALSE) +
  geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  geom_vline(xintercept = meanpH, linetype = 'dotted') +
  theme(legend.position='top') +
  guides(colour=guide_legend(ncol=3,bycol =TRUE,title.position = 'left')) +
  xlab('pH') + ylab(expression(paste('CO'[2]*' (mmol m'^{-2}*d^{-1}*')')))

## for the null model:
## for co2modnull
N <- 200
null.pdat <- with(droplevels(regvarf),
                  data.frame(`pH_surface` = rep(seq(7, 11, length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N))))
null.pred <- predict(co2modnull, newdata = null.pdat, type = "iterms") 
whichCols <- grep("pH", colnames(null.pred))

null.predse <- predict(co2modnull, newdata = null.pdat, type = "iterms", se.fit = TRUE)
null.predse <- as.data.frame(null.predse$se.fit)
whichColsse <- grep("pH", colnames(null.predse))

null.pdat <- cbind(null.pdat, Fitted = null.pred[, whichCols], Fittedse = null.predse[,whichColsse])
null.pdat <- with(null.pdat, transform(null.pdat, Fittedplus = Fitted + Fittedse))
null.pdat <- with(null.pdat, transform(null.pdat, Fittedminus = Fitted - Fittedse))

shiftnull <- attr(predict(co2modnull, newdata = null.pdat, type = "iterms"), "constant")
null.pdatnorm <- null.pdat
null.pdatnorm <- with(null.pdatnorm, transform(null.pdatnorm, Fitted = Fitted + shiftnull, 
                                               Fittedplus = Fittedplus + shiftnull, 
                                               Fittedminus = Fittedminus + shiftnull))

labdatco2 <- data.frame(x = 7.5, y = meanco2 - 18, label = "mean flux")

nullplot <- ggplot(null.pdatnorm, aes(x = pH_surface, y = Fitted)) +
  geom_line() +
  theme_bw(base_size = 15) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
          #  show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  #geom_point(data=regvars, aes(x=pH_surface, y=co2Flux)) +
  geom_vline(xintercept = meanpH, linetype="dotted") +
  ylab(expression(paste(CO[2]~"flux (mmol"~"C "*m^{-2}*"d"^{-1}*')'))) + xlab('pH')

# plot something for demonstrating GAMs
lmplot <- ggplot(regvars, aes(x = pH_surface, y = co2Flux)) +
  geom_point() +
  theme_bw(base_size = 15) +
  stat_smooth(method = "lm", linetype = "dotted", col='black') +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())

nullplotnodots <- ggplot(null.pdatnorm, aes(x = pH_surface, y = Fitted)) +
  geom_line(linetype='dotted') +
  theme_bw(base_size = 15) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_point(data=regvars, aes(x=pH_surface, y=co2Flux)) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
gamdemo <- plot_grid(lmplot, nullplotnodots)
ggsave("../docs/private/gamdemo.png", gamdemo, width=10, height=4, units = 'in')

## for egmod.red2
## oxygen
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

oxy.pdat <- with(droplevels(regvarf2),
                  data.frame(`Oxygen_ppm` = rep(seq(min(`Oxygen_ppm`, na.rm = TRUE),
                                                    max(`Oxygen_ppm`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
oxy.pdat <- merge(oxy.pdat, lakeXbar)
oxy.pred <- predict(egmod.red2, newdata = oxy.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("Oxy", colnames(oxy.pred$fit))
whichColsSE <- grep("Oxy", colnames(oxy.pred$se.fit))
oxy.pdat <- cbind(oxy.pdat, Fitted = oxy.pred$fit[, whichCols], 
                  se.Fitted = oxy.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
oxy.pdat <- with(oxy.pdat, transform(oxy.pdat, Fittedplus = Fitted + se.Fitted))
oxy.pdat <- with(oxy.pdat, transform(oxy.pdat, Fittedminus = Fitted - se.Fitted))

shiftoxy <- attr(predict(egmod.red2, newdata = oxy.pdat, type = "iterms"), "constant")
oxy.pdatnorm <- oxy.pdat
oxy.pdatnorm <- with(oxy.pdatnorm, transform(oxy.pdatnorm, Fitted = Fitted + shiftoxy, 
                                             Fittedplus = Fittedplus + shiftoxy, 
                                             Fittedminus = Fittedminus + shiftoxy))
labdatoxy <- data.frame(x = c(2.5, 11), y = c(meanpH + 0.03, 8.4), label = c("mean pH", 'mean oxygen'))

oxyquants <- quantile(regvarf2$Oxygen_ppm, c(.05,.95), na.rm = TRUE)

oxyplot <- ggplot(oxy.pdatnorm, aes(x = Oxygen_ppm, y = Fitted)) +
  annotate("rect", xmin=oxyquants[1], xmax=oxyquants[2], ymin=-Inf, ymax=Inf, alpha = .2) +
  geom_line() +
  theme_bw(base_size = 12) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanO2, linetype="dotted") +
  xlab(expression(paste(O[2]~"("*"mg"~L^{-1}*")"))) + ylab('pH')

## for GPP
N <- 200
varWant <- c("Oxygen_ppm", "TDN_ug_L", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

GPP.pdat <- with(droplevels(regvarf2),
                 data.frame(`GPP_h` = rep(seq(min(`GPP_h`, na.rm = TRUE),
                                                   max(`GPP_h`, na.rm = TRUE),
                                                   length = N),
                                               nlevels(Lake)),
                            Lake = rep(levels(Lake), each = N),
                            Year = rep(2004, prod(nlevels(Lake), N)),
                            dummy = rep(0, prod(nlevels(Lake), N))))
GPP.pdat <- merge(GPP.pdat, lakeXbar)
GPP.pred <- predict(egmod.red2, newdata = GPP.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("GPP", colnames(GPP.pred$fit))
whichColsSE <- grep("GPP", colnames(GPP.pred$se.fit))
GPP.pdat <- cbind(GPP.pdat, Fitted = GPP.pred$fit[, whichCols], 
                  se.Fitted = GPP.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
GPP.pdat <- with(GPP.pdat, transform(GPP.pdat, Fittedplus = Fitted + se.Fitted))
GPP.pdat <- with(GPP.pdat, transform(GPP.pdat, Fittedminus = Fitted - se.Fitted))

shiftGPP <- attr(predict(egmod.red2, newdata = GPP.pdat, type = "iterms"), "constant")
GPP.pdatnorm <- GPP.pdat
GPP.pdatnorm <- with(GPP.pdatnorm, transform(GPP.pdatnorm, Fitted = Fitted + shiftGPP, 
                                             Fittedplus = Fittedplus + shiftGPP, 
                                             Fittedminus = Fittedminus + shiftGPP))
labdatGPP <- data.frame(x = 0, y = meanpH + 0.01, label = "mean pH")

GPPplot <- ggplot(GPP.pdatnorm, aes(x = GPP_h, y = Fitted)) +
  geom_line() +
  theme_bw(base_size = 15) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanGPP, linetype='dotted')
  xlab(expression(paste('GPP ('~O[2]~h^{-1}*")"))) + ylab('pH')

## for TDN
N <- 200
varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

TDN.pdat <- with(droplevels(regvarf2),
                 data.frame(`TDN_ug_L` = rep(seq(min(`TDN_ug_L`, na.rm = TRUE),
                                              max(`TDN_ug_L`, na.rm = TRUE),
                                              length = N),
                                          nlevels(Lake)),
                            Lake = rep(levels(Lake), each = N),
                            Year = rep(2004, prod(nlevels(Lake), N)),
                            dummy = rep(0, prod(nlevels(Lake), N))))
TDN.pdat <- merge(TDN.pdat, lakeXbar)
TDN.pred <- predict(egmod.red2, newdata = TDN.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("TDN", colnames(TDN.pred$fit))
whichColsSE <- grep("TDN", colnames(TDN.pred$se.fit))
TDN.pdat <- cbind(TDN.pdat, Fitted = TDN.pred$fit[, whichCols], 
                  se.Fitted = TDN.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedplus = Fitted + se.Fitted))
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedminus = Fitted - se.Fitted))

shiftTDN <- attr(predict(egmod.red2, newdata = TDN.pdat, type = "iterms"), "constant")
TDN.pdatnorm <- TDN.pdat
TDN.pdatnorm <- with(TDN.pdatnorm, transform(TDN.pdatnorm, Fitted = Fitted + shiftTDN, 
                                             Fittedplus = Fittedplus + shiftTDN, 
                                             Fittedminus = Fittedminus + shiftTDN))
labdatN <- data.frame(x = 3000, y = meanpH + 0.03, label = "mean pH")

TDNquants <- quantile(regvarf2$TDN_ug_L, c(.05,.95), na.rm = TRUE)

TDNplot <- ggplot(TDN.pdatnorm, aes(x = TDN_ug_L, y = Fitted)) +
  annotate("rect", xmin=TDNquants[1], xmax=TDNquants[2], ymin=-Inf, ymax=Inf, alpha = .2) +
  geom_line() +
  theme_bw(base_size = 12) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                                     alpha = 0.25) +  
  #geom_text(data = labdatN, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE) +
  scale_x_log10(breaks = c(100,200,400,800,1600,3200,6400)) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanTDN, linetype = 'dotted') +
  xlab(expression(paste("TDN ("~mu*"g"~L^{-1}*")"))) + ylab('pH')

## for chl a
## Get site-specific data/predictions
N <- 1000
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

chla.pdat <- with(droplevels(regvarf2),
                  data.frame(`Chl_a_ug_L` = rep(seq(min(`Chl_a_ug_L`, na.rm = TRUE),
                                                    max(`Chl_a_ug_L`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
chla.pdat <- merge(chla.pdat, lakeXbar)
chla.pred <- predict(egmod.red2, newdata = chla.pdat, type = "terms")
whichCols <- grep("Chl", colnames(chla.pred))
chla.pdat <- cbind(chla.pdat, Fitted = rowSums(chla.pred[, whichCols]))

shiftchl <- attr(predict(egmod.red2, newdata = chla.pdat, type = "iterms"), "constant")
chl.pdatnorm <- chla.pdat
chl.pdatnorm <- with(chl.pdatnorm, transform(chl.pdatnorm, Fitted = Fitted + shiftchl))
labdatchl <- data.frame(x = 100, y = meanpH + 0.04, label = "mean pH")

chlquants <- quantile(regvarf2$Chl_a_ug_L, c(.05,.95), na.rm = TRUE)
chlaplot <- ggplot(chl.pdatnorm, aes(x = Chl_a_ug_L, y = Fitted, group = Lake, colour = 
                        ifelse(Lake == "WW", "Wascana", 
                               ifelse(Lake == "B", "Buffalo Pound", 
                                      "Katepwa, Pasqua, Diefenbaker,\nLast Mountain, Crooked")))) +
  annotate("rect", xmin=chlquants[1], xmax=chlquants[2], ymin=-Inf, ymax=Inf, alpha = .2) +
  geom_line() + 
  #geom_text(data = labdatchl, aes(label = label, x = x, y = y, size = 5),
  #          show.legend = FALSE, inherit.aes = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanchl, linetype='dotted') +
  papertheme + 
  scale_colour_brewer(name = "Lake", type = 'qual', palette = 'Dark2', direction=1) +
  theme(legend.position = 'top', legend.direction = "vertical", 
        axis.text.x = element_text(angle = 45)) +
  scale_x_log10(breaks = c(5,10,25,50,100,200,300)) +
  xlab(expression(paste("Chl"~italic(a)~"("~mu*"g"~"L"^{-1}*")"))) + 
  guides(colour=guide_legend(ncol=2,bycol =TRUE,title.position = 'left')) +
  ylab('pH')
## see http://stackoverflow.com/questions/11838278/
##    plot-with-conditional-colors-based-on-values-in-r
## make this bit into soi-pdo
N <- 100
simSOI <- c(-1, 0.3, 1.1)
SOIgroup <- factor(rep(c('-1', '0.3', '1.1'), each=100, times=7))  # 7*300
reptimes <- length(simSOI) #3
preddf <- data.frame(PDO = rep(seq(min(regvarf2$PDO, na.rm=TRUE),
                                      max(regvarf2$PDO, na.rm=TRUE),
                                      length = N), times=reptimes),
                     SOI = rep(simSOI, each = N)) #300
superN <- nrow(preddf)

varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

SOI.pdat <- with(droplevels(regvarf2),
                 data.frame(PDO = rep(preddf$PDO,
                                             nlevels(Lake)), # 300*7
                            SOI = rep(preddf$SOI,
                                        nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

SOI.pdat <- merge(SOI.pdat, lakeXbar)
SOI.pred <- predict(egmod.red2, newdata = SOI.pdat, type = "link")
SOI.pdat <- cbind(SOI.pdat, SOI.pred)
SOI.pdat$SOIgroup <- SOIgroup

names(SOI.pdat)[which(names(SOI.pdat)=='SOI.pred')] <- 'pH'

# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(SOI.pdat$PDO, SOI.pdat$SOI, regvarf2$PDO, regvarf2$SOI, dist=0.1)
SOI.pdat$pH[toofar] <- NA

labdatSOI <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
SOI.pdat$SOIgroup <- factor(SOI.pdat$SOIgroup, levels = c('-1', '0.3', '1.1'))

SOIplot <- ggplot(SOI.pdat[SOI.pdat$Lake=='B',], aes(x = PDO, y = pH, group= SOI, col=SOIgroup)) +
  papertheme +
  geom_line() +
  scale_color_discrete(name='SOI') +
  #scale_colour_brewer(name = "SOI", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  xlab('PDO') + ylab('pH')

## now slices for PDO
N <- 100
simPDO <- c(-1.5, 0, 1.6)
PDOgroup <- factor(rep(c('-1.5', '0', '1.6'), each=100, times=7))  # 7*300
reptimes <- length(simPDO) #3
preddf <- data.frame(SOI = rep(seq(min(regvarf2$SOI, na.rm=TRUE),
                                   max(regvarf2$SOI, na.rm=TRUE),
                                   length = N), times=reptimes),
                     PDO = rep(simPDO, each = N)) #300
superN <- nrow(preddf)

varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

PDO.pdat <- with(droplevels(regvarf2),
                 data.frame(SOI = rep(preddf$SOI,
                                      nlevels(Lake)), # 300*7
                            PDO = rep(preddf$PDO,
                                      nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

PDO.pdat <- merge(PDO.pdat, lakeXbar)
PDO.pred <- predict(egmod.red2, newdata = PDO.pdat, type = "link")
PDO.pdat <- cbind(PDO.pdat, PDO.pred)
PDO.pdat$PDOgroup <- PDOgroup

names(PDO.pdat)[which(names(PDO.pdat)=='PDO.pred')] <- 'pH'

# need to take away places where PDO*PDO combo has not occurred in the data!!
toofar <- exclude.too.far(PDO.pdat$PDO, PDO.pdat$SOI, regvarf2$PDO, regvarf2$SOI, dist=0.1)
PDO.pdat$pH[toofar] <- NA

labdatPDO <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
PDO.pdat$PDOgroup <- factor(PDO.pdat$PDOgroup, levels = c('-1.5', '0', '1.6'))

PDOplot <- ggplot(PDO.pdat[PDO.pdat$Lake=='B',], aes(x = SOI, y = pH, group= PDO, col=PDOgroup)) +
  papertheme +
  geom_line() +
  scale_color_discrete(name='PDO') +
  #scale_colour_brewer(name = "PDO", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH) +
  xlab('SOI') + ylab('pH')

## get expand grid object to plot a ggplot contour/heatmap
preddf <- with(regvarf2, expand.grid(PDO = seq(min(PDO), max(PDO), length = 100), 
                                   SOI = seq(min(SOI), max(SOI), length = 100)))
superN <- nrow(preddf)

varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

comb.pdat <- with(droplevels(regvarf2),
                 data.frame(PDO = rep(preddf$PDO,
                                        nlevels(Lake)), 
                            SOI = rep(preddf$SOI,
                                        nlevels(Lake)),
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

comb.pdat <- merge(comb.pdat, lakeXbar, sort=FALSE)
comb.pred <- predict(egmod.red2, newdata = comb.pdat, type = "terms")

whichCols <- grep("SOI", colnames(comb.pred))
comb.pdat <- cbind(comb.pdat, Fitted = comb.pred[, whichCols])

shiftcomb <- attr(comb.pred, "constant")
comb.pdatnorm <- comb.pdat
comb.pdatnorm <- with(comb.pdatnorm, transform(comb.pdatnorm, Fitted = Fitted + shiftcomb))

#exclude things too far from real
toofar <- exclude.too.far(comb.pdatnorm$PDO, comb.pdatnorm$SOI, regvarf2$PDO, regvarf2$SOI, dist=0.1)
comb.pdatnorm$pH <- comb.pdatnorm$Fitted
comb.pdatnorm$pH[toofar] <- NA

names(comb.pdat)[which(names(comb.pdat)=='SOI.pred')] <- 'pH'
comboplot <- ggplot(comb.pdatnorm, aes(x = SOI, y = PDO, z=pH)) + #, z=Fitted
  papertheme +
  geom_raster(aes(fill=pH)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent') +
  geom_point(data=regvarf2, aes(x=SOI, y=PDO, z=NULL)) +
  theme(legend.key.width=unit(2,"cm")) +
  geom_vline(linetype='dotted', xintercept = c(-1,0.3, 1.1)) +
  geom_abline(slope = 0, intercept = c(-1.5, 0, 1.6), linetype='dashed')

## save objects:
plots <- list(TDNplot, oxyplot, GPPplot, chlaplot, SOIplot, PDOplot, comboplot, nullplot, co2plot)
plotnames <- c('TDNplot', 'oxyplot', 'GPPplot','chlaplot', 'SOIplot', 'PDOplot', 'comboplot',
               'nullplot','co2plot')

invisible( # this means I don't get the list [[1:3]] returned on screen
  lapply(
    seq_along(plots), 
    function(x) ggsave(filename=paste0("../docs/private/gam-plots-", plotnames[x], ".pdf"), 
                       plot=plots[[x]], scale=0.8) # width=7, height=5, units = 'in'
  ) )

## arrange plots
## ## FIXME: get relative sizes of plots!!
allgam <- plot_grid(co2plot, chlaplot, oxyplot, TDNplot, ncol = 2, rel_heights = c(2,1.5), 
                    labels='AUTO')

climgam <- grid.arrange(comboplot, SOIplot, PDOplot, ncol = 2, layout_matrix = cbind(c(1,1,2), c(1,1,3)))

ggsave("../docs/private/ph-allgams.pdf", allgam, scale=0.77) #width=28, height=18, units = 'cm'
ggsave("../docs/private/climgam.pdf", climgam, scale=0.77, width = 7.5)
