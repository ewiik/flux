## These are the defaults for plotting my paper figures

## packages
library('ggplot2')
library('viridis')

## plot background: white
theme_bw() 
  
## viridis colour scheme
scale_fill_viridis()
scale_color_viridis()
# for discrete/factor variables
scale_color_viridis(discrete = TRUE, option = 'plasma')

## symbols and italics
expression(paste(italic('p')*'CO'[2]))
ylab(expression(paste(mu*'atm')))
~h^{-1}
expression(paste("VAR ("~mu*"g"~L^{-1}*")"))
expression(paste(CO[2]~"flux (mmol"~"C "*m^{-2}*"d"^{-1}*')'))

## ggplot labels for strips
varnames <- list(
  'Temperature'=expression(paste("Temperature ("~degree*"C)")) ,
  'Conductivity'= expression(paste("Conductivity ("*mu*"S"~"cm"^{-1}*")s")),
  'pH'="pH",
  'meanWindMS'= expression(paste("Mean Wind (m"~"s"^{-1}*")")),
  'SalCalc' = "Salinity (ppt)",
  'TICumol' = expression(paste("DIC ("~mu*"mol"~"L"^{-1}*")")),
  'Pressure' = 'Pressure (kPa)',
  'pco2atm' = expression(paste(italic(p)*"CO"[2]~"(ppm)"))
)

## how to get independent labeling to ggplot


  

