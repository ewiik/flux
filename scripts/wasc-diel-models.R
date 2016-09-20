## models for wascana data

## load packages
library("mgcv")
library("ggplot2")
library("viridis")

## load data
if (!file.exists('../data/private/WS-9-9-mod.rds')) {
  source("../scripts/dieltrial.R")}
wasc <- readRDS('../data/private/WS-9-9-mod.rds')

## run model for temporal patterns
m0 <- gam(pH ~ ti(Time, bs = "cc") + ti(Day), 
          data = wasc)
m1 <- gam(pH ~ ti(Time, bs = "cc") + ti(Day) + ti(Time, Day, bs = c("cc","tp")), 
          data = wasc)
anova(m0, m1, test = "LRT") # m1 better

m2 <- gam(pH ~ te(Time, Day, bs = c("cc","tp")), data=wasc)

## save models
saveRDS(m2, "../data/private/diel-wasc-timemod.rds")

## predict for periods of varying interactions.
N <- 200
simDOY <- c(254:256, 262:264, 270:272)
DOYgroup <- factor(rep(c('254:256','262:264','270:272'), each=3))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(wasc$Time, na.rm=TRUE),max(wasc$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
m2pred <- predict(m2, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, m2pred, DOYgroups)
names(predicted)[which(names(predicted)=='m2pred')] <- 'pH'

ggplot(predicted, aes(x = Time, y = pH, group= DOY, col=DOYgroups)) +
  papertheme +
  geom_line() +
  scale_colour_discrete(name="Day of Year") +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab('pH')


