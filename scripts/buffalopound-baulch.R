### quick look at Buffalo Pound monitoring data from 1985
## base file based on xlsx in fromkerri/

buff <- read.csv("data/private/BuffaloP_Data1985_2011.csv", na.strings = c("<NA>", "", "<"), 
                 as.is = "ODOsat")
buff$ODOsat
buffs <- buff
toreplace <- which(buff$ODOsat == NA | buff$ODOsat == "ERR")
buff$ODOsat[toreplace] <- "NA"
buff$ODOsat <- as.numeric(buff$ODOsat)
buff <- buff[-which(buff$ODOsat > 200),]

with(buff, plot(pH ~ ODOsat, col = week))

bpsal <- subset(params, Lake == "B")
     