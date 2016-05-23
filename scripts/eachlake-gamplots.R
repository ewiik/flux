### Get site specific data/predictions
### FIXME: load in objects properly
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
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
ph.pred <- predict(co2mod, newdata = ph.pdat, type = "terms")
whichCols <- grep("pH", colnames(ph.pred))
ph.pdat <- cbind(ph.pdat, Fitted = rowSums(ph.pred[, whichCols]))

ggplot(ph.pdat, aes(x = pH_surface, y = Fitted, 
                    colour = ifelse(Lake == "L", "L", 
                                    ifelse(Lake == "B","B",
                                           ifelse(Lake == "C", "C",
                                                  ifelse(Lake == "D", "D", "others")))))) +
  geom_line() +
  scale_colour_manual(name = 'Lakes', values = c("#9900CC","#339900", "#E69F00",
                                                 "#CC9999","#000000")) +
  xlab('pH') + ylab('mean-0 CO2 flux')

## for Oxygen
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


ggplot(oxy.pdat, aes(x = Oxygen_ppm, y = Fitted)) +
  geom_line() +
  geom_errorbar(limits) +
  xlab('Oxygen (mg/L)') + ylab('mean-0 pH')

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


ggplot(GPP.pdat, aes(x = GPP_h, y = Fitted)) +
  geom_line() +
  geom_errorbar(limits) +
  xlab('GPP/h') + ylab('mean-0 pH')

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


ggplot(TDN.pdat, aes(x = TDN_ug_L, y = Fitted)) +
  geom_line() +
  geom_errorbar(limits) +
  xlab('TDN (ug/L)') + ylab('mean-0 pH')

## FIXME: add to original units as per this:
## transform back into original pH value for clarity:
# get mean/intercept pH used by the model
shiftph <- attr(predict(phwwmod, newdata = ww.pdat, type = "iterms"), "constant")
ww.pdatnorm <- ww.pdat
ww.pdatnorm <- with(ww.pdatnorm, transform(ww.pdatnorm, Fitted = Fitted + shiftph, 
                                           Fittedplus = Fittedplus + shiftph, 
                                           Fittedminus = Fittedminus + shiftph))

labdat2 <- data.frame(x = 270, y = 9.1, label = "mean pH: 9.1")
# without this step, the resolution of the text is really off for some reason
