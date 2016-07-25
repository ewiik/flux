### Get site specific data/predictions

## load packages
library('ggplot2')
library('viridis')
library('cowplot')
library('mgcv')
library('extrafont')
library('gridExtra')
library("reshape2")
library('grid')

## Set defaults
theme_set(theme_bw())

## load in data
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")
if (!file.exists("../data/weathers.rds")) {
  source("../scripts/climate-weather-modeling.R")
}
weathers <- readRDS('../data/weathers.rds')

## Load in gam models
egmodlaggedsimp <- readRDS("../data/private/egmodlaggedsimp.rds")

source("../functions/geom_rug3.R")

## change to match what was used to create the models
regvars <- merge(regvars, weathers)
regvarf <- regvars
regvarf <- transform(regvarf, Year = as.factor(Year)) # make Year into factor for re
regvarf2 <- regvarf
regvarf2$`Chl_a_ug_L`[regvarf2$`Chl_a_ug_L` <=0 ] <- 0
regvarf2$`DOC_mg_L`[regvarf2$`DOC_mg_L` <=0 ] <- 0
regvarf2$`Chl_a_ug_L` <- regvarf2$`Chl_a_ug_L` + 1
regvarf2$`DOC_mg_L` <- regvarf2$`DOC_mg_L` + 1
regvarf2 <- transform(regvarf2, dummy = rep(1, nrow(regvarf2)))

## what is mean pH? perhaps add a dotted line at this point in the plots
meanpH <- with(regvarf2, 
               mean(pH_surface, na.rm=TRUE)) 
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
meanspei <- with(regvarf2, 
                 mean(SPEI02, na.rm=TRUE))
## limit prediction data frame to observed intervals
regsplit <- with(regvarf2, split(regvarf2, list(Lake)))
regsplit <- regsplit[sapply(regsplit, function(x) dim(x)[1]) > 0] #remove empties for lakes R, E etc.
minmax <- function(df, colnames) {
  allmin <- as.data.frame(do.call(cbind, lapply(df[,colnames], min, na.rm = TRUE)))
  names(allmin) <- sapply(names(allmin), function(x) paste("min",x, sep = ""))
  allmax <- as.data.frame(do.call(cbind, lapply(df[,colnames], max, na.rm = TRUE)))
  names(allmax) <- sapply(names(allmax), function(x) paste("max",x, sep = ""))
  summ <- as.data.frame(cbind(allmin, allmax))
  summ <- cbind(summ, data.frame(Lake = df$Lake[1]))
  summ
}
minmaxes <- do.call(rbind, lapply(regsplit, minmax, colnames = c("GPP_h", "TDN_ug_L", "DOC_mg_L", 
                                                                 "Oxygen_ppm", "pH_surface", "Chl_a_ug_L")))
rownames(minmaxes) <- NULL

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=18, base_family = 'Arial') +
  theme(legend.position='top')


## for SPEI
N <- 200
varWant <- c("TDN_ug_L", "Chl_a_ug_L", "PDOmean", "SOImean")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

spei.pdat <- with(droplevels(regvarf2),
                  data.frame(`SPEI02` = rep(seq(min(`SPEI02`, na.rm = TRUE),
                                                max(`SPEI02`, na.rm = TRUE),
                                                length = N),
                                            nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
spei.pdat <- merge(spei.pdat, lakeXbar)
spei.pred <- predict(egmodlaggedsimp, newdata = spei.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("SPEI", colnames(spei.pred$fit))
whichColsSE <- grep("SPEI", colnames(spei.pred$se.fit))
spei.pdat <- cbind(spei.pdat, Fitted = spei.pred$fit[, whichCols], 
                   se.Fitted = spei.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
spei.pdat <- with(spei.pdat, transform(spei.pdat, Fittedplus = Fitted + se.Fitted))
spei.pdat <- with(spei.pdat, transform(spei.pdat, Fittedminus = Fitted - se.Fitted))

shiftspei <- attr(predict(egmodlaggedsimp, newdata = spei.pdat, type = "iterms"), "constant")
spei.pdatnorm <- spei.pdat
spei.pdatnorm <- with(spei.pdatnorm, transform(spei.pdatnorm, Fitted = Fitted + shiftspei, 
                                               Fittedplus = Fittedplus + shiftspei, 
                                               Fittedminus = Fittedminus + shiftspei))

speiquants <- quantile(regvarf2$SPEI02, c(.05,.95), na.rm = TRUE)

speiplot <- ggplot(spei.pdatnorm, aes(x = SPEI02, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=speiquants[1], xmax=speiquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanspei, linetype="dotted") +
  xlab("SPEI index (dry to wet)") + ylab('pH') +
  ylim(c(8.7, 9.2))

## for TDN
N <- 200
varWant <- c("Chl_a_ug_L", "PDOmean", "SOImean", "SPEI02")
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
TDN.pred <- predict(egmodlaggedsimp, newdata = TDN.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("TDN", colnames(TDN.pred$fit))
whichColsSE <- grep("TDN", colnames(TDN.pred$se.fit))
TDN.pdat <- cbind(TDN.pdat, Fitted = TDN.pred$fit[, whichCols], 
                  se.Fitted = TDN.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedplus = Fitted + se.Fitted))
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedminus = Fitted - se.Fitted))

shiftTDN <- attr(predict(egmodlaggedsimp, newdata = TDN.pdat, type = "iterms"), "constant")
TDN.pdatnorm <- TDN.pdat
TDN.pdatnorm <- with(TDN.pdatnorm, transform(TDN.pdatnorm, Fitted = Fitted + shiftTDN, 
                                             Fittedplus = Fittedplus + shiftTDN, 
                                             Fittedminus = Fittedminus + shiftTDN))
TDN.pdatnorm <- merge(TDN.pdatnorm, minmaxes)
overs <- with(TDN.pdatnorm, which(TDN_ug_L < minTDN_ug_L | TDN_ug_L > maxTDN_ug_L))
TDN.pdatnorm <- TDN.pdatnorm[-overs,]

labdatN <- data.frame(x = 3000, y = meanpH + 0.03, label = "mean pH")

TDNquants <- quantile(regvarf2$TDN_ug_L, c(.05,.95), na.rm = TRUE)

TDNplot <- ggplot(TDN.pdatnorm, aes(x = TDN_ug_L, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=TDNquants[1], xmax=TDNquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_rug3(aes(x=TDN_ug_L), data = regvarf2, stat = "identity", position = "identity", 
            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  #geom_text(data = labdatN, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  scale_x_log10(breaks = c(100,200,400,800,1600,3200,6400)) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanTDN, linetype = 'dotted') +
  xlab(expression(paste("TDN ("~mu*"g"~L^{-1}*")"))) + ylab('pH')

## for chl a
## Get site-specific data/predictions
N <- 1000
varWant <- c("TDN_ug_L", "PDOmean", "SOImean", "SPEI02")
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
chla.pred <- predict(egmodlaggedsimp, newdata = chla.pdat, type = "terms")
whichCols <- grep("Chl", colnames(chla.pred))
chla.pdat <- cbind(chla.pdat, Fitted = rowSums(chla.pred[, whichCols]))

shiftchl <- attr(predict(egmodlaggedsimp, newdata = chla.pdat, type = "iterms"), "constant")
chl.pdatnorm <- chla.pdat
chl.pdatnorm <- with(chl.pdatnorm, transform(chl.pdatnorm, Fitted = Fitted + shiftchl))
chl.pdatnorm <- merge(chl.pdatnorm, minmaxes)
overs <- with(chl.pdatnorm, which(Chl_a_ug_L < minChl_a_ug_L | Chl_a_ug_L > maxChl_a_ug_L))
chl.pdatnorm <- chl.pdatnorm[-overs,]

labdatchl <- data.frame(x = 100, y = meanpH + 0.04, label = "mean pH")

chlquants <- quantile(regvarf2$Chl_a_ug_L, c(.05,.95), na.rm = TRUE)
chlaplot <- ggplot(chl.pdatnorm, aes(x = Chl_a_ug_L, y = Fitted, group = Lake, colour = 
                                       ifelse(Lake == "WW", "Wascana", 
                                              ifelse(Lake == "D", "Diefenbaker",
                                                     ifelse(Lake == 'K', "Katepwa", 
                                                            ifelse(Lake == 'P', "Pasqua", 
                                                                   ifelse(Lake == 'B', 'Buffalo Pound', 
                                                                          ifelse(Lake=='L','Last Mountain',
                                                                                 'Crooked')))))),
                                     lty=ifelse(Lake == "WW", "Wascana", 
                                                ifelse(Lake == "D", "Diefenbaker",
                                                       ifelse(Lake == 'K', "Katepwa", 
                                                              ifelse(Lake == 'P', "Pasqua", 
                                                                     ifelse(Lake == 'B', 'Buffalo Pound', 
                                                                            ifelse(Lake=='L','Last Mountain',
                                                                                   'Crooked')))))))) +
  papertheme + 
  annotate("rect", xmin=chlquants[1], xmax=chlquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() + 
  #geom_text(data = labdatchl, aes(label = label, x = x, y = y, size = 5),
  #          show.legend = FALSE, inherit.aes = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanchl, linetype='dotted') +
  scale_linetype_manual(name='Lake', values = c("solid", "solid","dotdash", "solid", "solid", 
                                                "longdash", "longdash")) +
  scale_colour_manual(name="Lake", values = c("#5e3c99", "#5e3c99","#5e3c99", "#b2abd2", "#5e3c99",
                                              "#b2abd2", "#e66101"))+
  theme(legend.position = 'top', legend.direction = "vertical", 
        axis.text.x = element_text(angle = 45)) +
  scale_x_log10(breaks = c(5,10,25,50,100,200,300)) +
  xlab(expression(paste("Chl"~italic(a)~"("~mu*"g"~"L"^{-1}*")"))) + 
  guides(colour=guide_legend(ncol=4,nrow=2,bycol =TRUE,title.position = 'left'),
         lty=guide_legend(ncol=4,nrow=2, bycol =TRUE,title.position = 'left')) +
  ylab('pH')
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

## see http://stackoverflow.com/questions/11838278/
##    plot-with-conditional-colors-based-on-values-in-r
## make this bit into soi-pdo
N <- 100
simSOI <- c(-1.1, 0.3, 1.1)
SOIgroup <- factor(rep(c('-1.1', '0.3', '1.1'), each=100, times=7))  # 7*300
reptimes <- length(simSOI) #3
preddf <- data.frame(PDOmean = rep(seq(min(regvarf2$PDOmean, na.rm=TRUE),
                                   max(regvarf2$PDOmean, na.rm=TRUE),
                                   length = N), times=reptimes),
                     SOImean = rep(simSOI, each = N)) #300
superN <- nrow(preddf)

varWant <- c("SPEI02", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

SOI.pdat <- with(droplevels(regvarf2),
                 data.frame(PDOmean = rep(preddf$PDOmean,
                                      nlevels(Lake)), # 300*7
                            SOImean = rep(preddf$SOImean,
                                      nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

SOI.pdat <- merge(SOI.pdat, lakeXbar)
SOI.pred <- predict(egmodlaggedsimp, newdata = SOI.pdat, type = "link")
SOI.pdat <- cbind(SOI.pdat, SOI.pred)
SOI.pdat$SOIgroup <- SOIgroup

names(SOI.pdat)[which(names(SOI.pdat)=='SOI.pred')] <- 'pH'

# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(SOI.pdat$PDOmean, SOI.pdat$SOImean, 
                          regvarf2$PDOmean, regvarf2$SOImean, dist=0.1)
SOI.pdat$pH[toofar] <- NA

labdatSOI <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
SOI.pdat$SOIgroup <- factor(SOI.pdat$SOIgroup, levels = c('-1.1', '0.3', '1.1'))

SOIplot <- ggplot(SOI.pdat[SOI.pdat$Lake=='B',], aes(x = PDOmean, y = pH, group= SOImean, col=SOIgroup, 
                                                     lty=SOIgroup)) +
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_line() +
  scale_color_manual(name=expression(paste(bold('E')~'      SOI')), values = c("#5e3c99", "#b2abd2", "#e66101"))+
  scale_linetype_manual(name=expression(paste(bold('E')~'      SOI')), values = c("solid", "solid","longdash")) +
  #scale_colour_brewer(name = "SOI", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  xlab('PDO') + ylab('pH')
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

## now slices for PDO
N <- 100
simPDO <- c(-0.5, 0.5, 1.5)
PDOgroup <- factor(rep(c('-0.5', '0.5', '1.5'), each=100, times=7))  # 7*300
reptimes <- length(simPDO) #3
preddf <- data.frame(SOImean = rep(seq(min(regvarf2$SOImean, na.rm=TRUE),
                                   max(regvarf2$SOImean, na.rm=TRUE),
                                   length = N), times=reptimes),
                     PDOmean = rep(simPDO, each = N)) #300
superN <- nrow(preddf)

varWant <- c("SPEI02", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

PDO.pdat <- with(droplevels(regvarf2),
                 data.frame(SOImean = rep(preddf$SOImean,
                                      nlevels(Lake)), # 300*7
                            PDOmean = rep(preddf$PDOmean,
                                      nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

PDO.pdat <- merge(PDO.pdat, lakeXbar)
PDO.pred <- predict(egmodlaggedsimp, newdata = PDO.pdat, type = "link")
PDO.pdat <- cbind(PDO.pdat, PDO.pred)
PDO.pdat$PDOgroup <- PDOgroup

names(PDO.pdat)[which(names(PDO.pdat)=='PDO.pred')] <- 'pH'

# need to take away places where PDO*PDO combo has not occurred in the data!!
toofar <- exclude.too.far(PDO.pdat$PDOmean, PDO.pdat$SOImean, regvarf2$PDOmean, regvarf2$SOImean, dist=0.1)
PDO.pdat$pH[toofar] <- NA

labdatPDO <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
PDO.pdat$PDOgroup <- factor(PDO.pdat$PDOgroup, levels = c('-0.5', '0.5', '1.5'))

PDOplot <- ggplot(PDO.pdat[PDO.pdat$Lake=='B',], aes(x = SOImean, y = pH, group= PDOmean, col=PDOgroup,
                                                     lty=PDOgroup)) +
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_line() +
  scale_color_manual(name=expression(paste(bold('F')~'      PDO')), values = c("#5e3c99", "#b2abd2", "#e66101"))+
  scale_linetype_manual(name=expression(paste(bold('F')~ '      PDO')), values = c("solid", "solid","longdash")) +
  #scale_colour_brewer(name = "PDO", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH) +
  xlab('SOI') + ylab('pH')

## get expand grid object to plot a ggplot contour/heatmap
preddf <- with(regvarf2, expand.grid(PDOmean = seq(min(PDOmean), max(PDOmean), length = 100), 
                                     SOImean = seq(min(SOImean), max(SOImean), length = 100)))
superN <- nrow(preddf)

varWant <- c("SPEI02", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

comb.pdat <- with(droplevels(regvarf2),
                  data.frame(PDOmean = rep(preddf$PDOmean,
                                       nlevels(Lake)), 
                             SOImean = rep(preddf$SOImean,
                                       nlevels(Lake)),
                             Lake = rep(levels(Lake), each = superN),
                             Year = rep(2004, prod(nlevels(Lake), superN)),
                             dummy = rep(0, prod(nlevels(Lake), superN))))

comb.pdat <- merge(comb.pdat, lakeXbar, sort=FALSE)
comb.pred <- predict(egmodlaggedsimp, newdata = comb.pdat, type = "terms")

whichCols <- grep("SOI", colnames(comb.pred))
comb.pdat <- cbind(comb.pdat, Fitted = comb.pred[, whichCols])

shiftcomb <- attr(comb.pred, "constant")
comb.pdatnorm <- comb.pdat
comb.pdatnorm <- with(comb.pdatnorm, transform(comb.pdatnorm, Fitted = Fitted + shiftcomb))

#exclude things too far from real
toofar <- exclude.too.far(comb.pdatnorm$PDOmean, comb.pdatnorm$SOImean, regvarf2$PDOmean, regvarf2$SOImean, dist=0.1)
comb.pdatnorm$pH <- comb.pdatnorm$Fitted
comb.pdatnorm$pH[toofar] <- NA

names(comb.pdat)[which(names(comb.pdat)=='SOI.pred')] <- 'pH'
comboplot <- ggplot(comb.pdatnorm, aes(x = SOImean, y = PDOmean, z=pH)) + #, z=Fitted
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=pH)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent') +
  geom_point(data=regvarf2, aes(x=SOImean, y=PDOmean, z=NULL)) +
  theme(legend.key.width=unit(2,"cm")) +
  geom_vline(linetype='dashed', xintercept = c(-1.1, 0.3, 1.1)) +
  geom_abline(slope = 0, intercept = c(-0.5, 0.5, 1.5), linetype='dashed') +
  ylab("PDO") + xlab("SOI")

## save objects:
plots <- list(TDNplot, chlaplot, SOIplot, PDOplot, comboplot, speiplot)
plotnames <- c('TDNplot', 'chlaplot', 'SOIplot', 'PDOplot', 'comboplot', 'SPEIplot')

invisible( # this means I don't get the list [[1:3]] returned on screen
  lapply(
    seq_along(plots), 
    function(x) ggsave(filename=paste0("../docs/private/gam-plots-lagged-", plotnames[x], ".pdf"), 
                       plot=plots[[x]], scale=0.6) # width=7, height=5, units = 'in'
  ) )

## arrange plots
allgam <- plot_grid(chlaplot, TDNplot, speiplot, ncol = 1, rel_heights = c(1.5,1,1),  labels='AUTO')
climgam <- grid.arrange(comboplot, SOIplot, PDOplot, ncol = 2, layout_matrix = cbind(c(1,1,2), c(1,1,3)))
#papergam <- grid.arrange(climgam, chlaplot, oxyplot, TDNplot, 
#                     layout_matrix=rbind(c(2,1),c(2,1), c(2,1), c(3,1), c(3,1), c(4,1), c(4,1)))
## FIXME: how to get labels in grid.arrange??
splineplots <- plot_grid(chlaplot, TDNplot, speiplot, ncol = 1, nrow=3, rel_heights = c(2, 1.5, 1.5), 
                         labels='AUTO')
vararrange <- plot_grid(splineplots, climgam, ncol = 2, labels=c('', "D"))
ggsave("../docs/private/ph-allgams-lagged.pdf", allgam, scale=0.77) #width=28, height=18, units = 'cm'
ggsave("../docs/private/climgam-lagged.pdf", climgam, scale=0.77, width = 7.5)
ggsave("../docs/private/paper-gams-lagged.pdf", vararrange, scale=1.2, width = 12)

