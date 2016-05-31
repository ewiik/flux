## plotting diel data together

## read in data
ww <- readRDS('../data/private/WS-9-9-mod.rds')
bp <- readRDS('../data/private/bpbuoy2014-mod.rds')

## add week to data
ww <- transform(ww, Week = as.numeric(format(Date.Time, "%V")))

## subset to what we wanna merge
wwsub <- subset(ww, select = c('pCO2', 'isDay'))
bpsub <- subset(bp, select=c('co2corr', 'isDay'))

## add lake designator and change colnames
wwsub$Lake <- rep('Wascana')
bpsub$Lake <- rep('BuffaloPound')

names(bpsub)[which(names(bpsub) == 'co2corr')] <- 'pCO2'

## rbind them
alldat <- rbind(wwsub, bpsub)

## change isDay to something more plottable
nights <- which(alldat$isDay == FALSE)
days <- which(alldat$isDay == TRUE)
  
alldat$isDay[nights] <- 'Night'
alldat$isDay[days] <- 'Day'
alldat$isDay <- as.factor(alldat$isDay)
isDay <- alldat$isDay
alldat$isDay <- factor(alldat$isDay,levels(alldat$isDay)[c(2,1)])


## plot
ggplot(alldat, aes(isDay, pCO2, fill=Lake)) +
  geom_boxplot() +
  geom_hline(yintercept = c(2000,400)) +
  geom_text(x=1.5, y=450, label='~atmospheric pCO2') +
  geom_text(x=1.5, y=1950, label='probe limit') +
  xlab(label = '')

ggplot(bp, aes(factor(Hour), co2corr, col=factor(Month), fill=factor(Month))) +
  geom_boxplot(outlier.colour = NULL) +
  geom_hline(yintercept = c(2000,400)) +
  #geom_text(x=1.5, y=450, label='~atmospheric pCO2') +
  geom_text(x=5, y=1950, label='probe limit', col="black") +
  xlab(label = 'Hour')

ggplot(ww, aes(factor(Hour), pCO2, col=factor(Week), fill=factor(Week))) +
  geom_boxplot(outlier.colour = NULL) +
  geom_hline(yintercept = 400) +
  #geom_text(x=1.5, y=450, label='~atmospheric pCO2') +
  #geom_text(x=5, y=1950, label='probe limit', col="black") +
  xlab(label = 'Hour')

## signif test
wilcox.test(bpsub$pCO2, wwsub$pCO2, alternative = 'two.sided')
ansari.test(bpsub$pCO2, wwsub$pCO2, alternative = 'two.sided')

wilcox.test(bpsub$pCO2[bpsub$isDay == TRUE], 
            bpsub$pCO2[bpsub$isDay == FALSE], alternative = 'two.sided')
wilcox.test(wwsub$pCO2[wwsub$isDay == TRUE], 
            wwsub$pCO2[wwsub$isDay == FALSE], alternative = 'two.sided')




