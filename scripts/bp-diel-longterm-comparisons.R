## compare diel and longterm for bp

## read in data
blong <- readRDS('../data/private/bp-longterm-mod.rds')
bdiel <-  readRDS('../data/private/bpbuoy2014-mod.rds')

blong <- subset(blong, year == 2014)
isna <- which(is.na(blong$pCO2uatm))

blong <- blong[-isna,]
