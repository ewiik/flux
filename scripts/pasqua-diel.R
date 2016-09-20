## pasqua diel; from Peter's email 29.8.2016
## 'bottom' = 8m; "top" = 1m; PAR depth pends checking by Jay
## units: µS/cm for cond, mg/L for ODO, µmol/s/m² for PAR

## load packages
library("ggplot2")
library("reshape2")
library("gridExtra")

## read in data and transform to correct format
pas <- read.csv("../data/private/PasquaAug25_2016Buoy.csv", skip=4, header = FALSE)
names(pas) <- c('datetime', 'battery', 'toptemp', 'topcond','topodosat','topodo',
                'bottomtemp','bottomcond','bottomodosat','bottomodo','PARtop')
pas <- transform(pas, datetime=as.POSIXct(as.character(datetime), format="%m/%d/%Y %H:%M"))

## melt for various purposes
oxymelt <- melt(pas[,c("datetime", "topodosat", "bottomodosat")], id.vars = "datetime")
oxymelt$variable <- as.character(oxymelt$variable)
oxymelt$variable[grep("top", oxymelt$variable)] <- "1m"
oxymelt$variable[grep("bott", oxymelt$variable)] <- "8m"
names(oxymelt)[grep("var", names(oxymelt))] <- "Depth"

tempmelt <- melt(pas[,c("datetime", "toptemp", "bottomtemp")], id.vars = "datetime")
tempmelt$variable <- as.character(tempmelt$variable)
tempmelt$variable[grep("top", tempmelt$variable)] <- "1m"
tempmelt$variable[grep("bott", tempmelt$variable)] <- "8m"
names(tempmelt)[grep("var", names(tempmelt))] <- "Depth"

## plot 
oxyplot <- ggplot(data=oxymelt, mapping = aes(y=value, x=datetime, group=Depth, lty=Depth)) +
  geom_line() +
  ylab("Oxygen saturation (%)") +
  xlab("Date")

tempplot <- ggplot(data=tempmelt, mapping = aes(y=value, x=datetime, group=Depth, lty=Depth)) +
  geom_line() +
  ylab(expression(paste("Temperature ("*degree*"C)"))) +
  xlab("Date")

comboplot <- grid.arrange(tempplot, oxyplot, ncol=1)
