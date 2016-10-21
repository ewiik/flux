## formulation of model for Gavin

gam(co2corr ~ s(chlrfu, by = factor(TimeofDay)) + s(bga1rfu, by = factor(TimeofDay)) + 
  s(airtemp) + s(stability) + s(ODOrel1) + s(dailyrain) + s(conv) + 
  s(turb) + s(cdom) + s(windsp), data=bdat, ...)