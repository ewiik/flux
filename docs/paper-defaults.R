## These are the defaults for plotting my paper figures
## assuming L&O guidelines:
## When possible, submit figures in the size you wish to have them appear in the journal. 
## Most illustrations, except some maps and very wide graphs, should be 1-column size (3.5 inches) 
## and a resolution of 300 dpi. The font size on the x and y axes should not be larger than that 
## of the title, and the same font (Arial or Times New Roman is preferred) should be used throughout. 
## Numbers on the x and y axes should be smaller than the descriptive title, which should be 
## 12-point font. Fonts smaller than 12 points are generally not legible when reduced to 
## 1 column size. Use boldface type with care; if illustrations are to be reduced, the letters 
## with open spaces will disappear. Use sentence case (capitalize the first word ONLY) 
## for axis titles, labels, and legends.
## http://aslopubs.onlinelibrary.wiley.com/hub/journal/10.1002/(ISSN)1939-5590/about/author-guidelines.html

## http://www.fromthebottomoftheheap.net/2013/09/09/preparing-figures-for-plos-one-with-r/

## packages
library('ggplot2')
library('viridis')
library('cowplot')
library('extrafont')

library('scales')
## plot background: white
theme_bw() 
  
## viridis colour scheme for contour plots
scale_fill_viridis()
scale_color_viridis()
# for discrete/factor variables
scale_color_viridis(discrete = TRUE, option = 'plasma')

## colour brewer scheme for regular plots
scale_colour_brewer(type = 'qual', ...)

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

## how get correct figure sizes and text
ggsave(.., width=3.5, units='in') # one column
# optional way to control scale of text vs plot: ggsave(.., scale=)
theme_bw(base_size = 12, base_family = 'Arial') # set baseline font size  
# optional level of control: only applies to legends and axis labels, not tick labels
# (..., text = element_text(size = ))

  
## play with colour selections...
show_col(brewer_pal("seq", palette='Oranges')(5))

