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
