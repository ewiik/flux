### Get site specific data/predictions
### FIXME: the co2 model still needs proper plotting, code needs adjusting for multiple lakes
###     scenario

## load in data
regvars <- readRDS("../data/private/regvars.rds")

## Load in gam models
co2mod <- readRDS("../data/private/co2mod.rds")
egmod.red2 <- readRDS("../data/private/egmodred2.rds")

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

## generate predicted for co2 and ph for all signif variables
## co2 ~ ph
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "phgen_ppm", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

ph.pdat <- with(droplevels(regvarf2),
                  data.frame(`pH_surface` = rep(seq(7,
                                                    11,
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
ph.pdat <- merge(ph.pdat, lakeXbar)
ph.pred <- predict(co2mod, newdata = ph.pdat, type = "terms", se.fit = TRUE)

whichCols <- grep("pH", colnames(ph.pred$fit))
whichColsSE <- grep("pH", colnames(ph.pred$se.fit))
ph.pdat <- cbind(ph.pdat, Fitted = ph.pred$fit[, whichCols], 
                  se.Fitted = ph.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
ph.pdat <- with(ph.pdat, transform(ph.pdat, Fittedplus = Fitted + se.Fitted))
ph.pdat <- with(ph.pdat, transform(ph.pdat, Fittedminus = Fitted - se.Fitted))

shiftph <- attr(predict(co2mod, newdata = ph.pdat, type = "iterms"), "constant")
ph.pdatnorm <- ph.pdat
ph.pdatnorm <- with(ph.pdatnorm, transform(ph.pdatnorm, Fitted = Fitted + shiftph, 
                                             Fittedplus = Fittedplus + shiftph, 
                                             Fittedminus = Fittedminus + shiftph))
labdatph <- data.frame(x = 2.5, y = meanpH + 0.03, label = "mean pH")

ggplot(ph.pdatnorm, aes(x = phgen_ppm, y = Fitted)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_text(data = labdatph, aes(label = label, x = x, y = y, size = 5), 
            show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  xlab('phgen (mg/L)') + ylab('pH')



ggplot(ph.pdat, aes(x = pH_surface, y = Fitted, 
                    colour = ifelse(Lake == "L", "L", 
                                    ifelse(Lake == "B","B",
                                           ifelse(Lake == "C", "C",
                                                  ifelse(Lake == "D", "D", "others")))))) +
  geom_line() +
  scale_colour_manual(name = 'Lakes', values = c("#9900CC","#339900", "#E69F00",
                                                 "#CC9999","#000000")) +
  xlab('pH') + ylab('mean-0 CO2 flux')

## for oxygen
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

oxy.pdat <- with(droplevels(regvarf2),
                  data.frame(`oxygen_ppm` = rep(seq(min(`oxygen_ppm`, na.rm = TRUE),
                                                    max(`oxygen_ppm`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
oxy.pdat <- merge(oxy.pdat, lakeXbar)
oxy.pred <- predict(egmod.red2, newdata = oxy.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("oxy", colnames(oxy.pred$fit))
whichColsSE <- grep("oxy", colnames(oxy.pred$se.fit))
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
labdatoxy <- data.frame(x = 2.5, y = meanoxy + 0.03, label = "mean oxy")

ggplot(oxy.pdatnorm, aes(x = Oxygen_ppm, y = Fitted)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_text(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
            show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  xlab('Oxygen (mg/L)') + ylab('pH')

## for GPP
N <- 200
varWant <- c("phgen_ppm", "TDN_ug_L", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
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

ggplot(GPP.pdatnorm, aes(x = GPP_h, y = Fitted)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_text(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
            show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  xlab('GPP/h') + ylab('pH')

## for TDN
N <- 200
varWant <- c("phgen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
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

ggplot(TDN.pdatnorm, aes(x = TDN_ug_L, y = Fitted)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                                     alpha = 0.25) +  
  geom_text(data = labdatN, aes(label = label, x = x, y = y, size = 5), 
            show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
    xlab('TDN (ug/L)') + ylab('pH')


