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

  

