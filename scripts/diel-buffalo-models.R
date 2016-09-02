## statistics on bp data

## read in data
if (!file.exists('../data/private/bpbuoy2014-mod.rds')) {
  source("diel-buffalo.R")}
bdat <- readRDS('../data/private/bpbuoy2014-mod.rds')

## stats
lm2014night <- with(bdat2014full[bdat2014full$isDay == FALSE,], lm(co2corr~pco2))
lm2014day <- with(bdat2014full[bdat2014full$isDay == TRUE,], lm(co2corr~pco2))
lm2014int <- with(bdat2014full, lm(co2corr~pco2))
with(bdat2014full, cor(x=co2corr, y=pco2, use='complete.obs', method='pearson'))

# just with the 10am-3pm data
ten3 <- subset(bdat2014full, Hour >= 10 & Hour <=15)

## could do this to make converge but also did manually, required four iterations
#reps <- 8 #takes a while for it to converge
#for (i in 1:reps) {
#   start <- if (i == 1 ) {coef(glm(co2corr ~ pco2, data = ten3, family = gaussian))}
#     else { coef(glmmod) }
#   glmmod <- glm(co2corr ~ pco2, data = ten3, family = Gamma(link = "identity"), start = start,
#                 control = glm.control(maxit=100))
# }

start <- c(-69.0709700, 0.9070629 )
glmmod <- glm(co2corr ~ pco2, data = ten3, family = Gamma(link = "identity"), start = start,
              control = glm.control(maxit=100))

pglm <- ggplot(data=ten3, aes(y = co2corr, x = pco2)) + 
  theme_bw() +
  geom_point() +
  geom_abline(slope=1, intercept=0, color='red', lty=3) +
  #geom_ribbon() + ## FIXME: add the regression line here
  geom_abline(intercept = glmmod$coefficients[1], slope=glmmod$coefficients[2], col = 'black') +
  scale_color_viridis(discrete = TRUE, 
                      option = 'plasma', direction = -1) +
  ylab(expression(paste('Measured'~italic('p')*'CO'[2]~'('*mu*'atm)'))) +
  xlab(expression(paste('Calculated'~italic('p')*'CO'[2]~'('*mu*'atm)')))


# residual plot etc
plot(resid(lm2014int) ~ lm2014int$fitted.values, ylab='Residuals (pCO2 uatm)', 
     xlab='Fitted values (pCO2 uatm)')

with(bdat2014full, plot(co2corr ~ pco2, col=ifelse(isDay, 'red', 'black'),
                        xlab = 'Calculated pCO2', ylab = 'Measured pCO2'))
abline(0,1, lty=3)
abline(lm2014day$coefficients[1], lm2014day$coefficients[2], col = 'red')
abline(lm2014night$coefficients[1], lm2014night$coefficients[2], col = 'black')
#abline(-71.6, 1.143408, lty = 3, col = 'green')
legend('topleft', col = c('red', 'black', 'black'), pch=c(1,1, NA), lty = c(1,1,3), 
       legend = c('daytime values and regression line', 
                  'night-time values and regression line',
                  '1:1 line'))
legend('bottomright', col = c('red', 'black'), title = 'R2', 
       legend = c('day = 0.72', 'night = 0.92'))

with(bdat2014full, plot(pco2 - co2corr ~ co2corr, xlab = 'Measured pCO2', 
                        ylab = 'Calculated - measured pCO2'))
abline(0,0)

## any DIC concentration evidence for carbonate precipitation?
ggplot(data=params[params$Lake =='B',], aes(y = TIC, x = Month, group = Year)) + # ifelse(Hour >= 8 & Hour <=20, 
  #   'red', 'black')
  #scale_color_identity() + # means it understands red and black in ifelse
  #scale_color_manual(values=c('black', 'red')) +
  geom_point() +
  facet_wrap('Year') +
  ylab("DIC (mg/L)") +
  xlab("Month")

## gamming the cyclic stuff in the data

m1 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(Week), data = bdat2014full)

m2 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(Week) + ti(Time, Week, bs = c("cc","tp")), 
          data = bdat2014full)
anova(m1, m2, test = "LRT")

m3 <- gam(ph1 ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")), 
          data = bdat2014full)
mco2 <- gam(co2corr ~ ti(Time, bs = "cc") + ti(DOY) + ti(Time, DOY, bs = c("cc","tp")), 
            data = bdat2014full)

plot(m2, scheme=2, ylab='mean-0 pH', main='mean-0 pH') # latter works for 3D plot
#   former for the 2D plot
plot(acf(resid(m2)))

pdf("diel-bpgam.pdf")
op <- par(mar = c(4,4,1,1) + 0.1)
plot(m2, pages = 4, scheme = 2)
par(op)
dev.off()

## predict for periods of varying interactions.
N <- 200
simDOY <- c(171:173, 211:213, 237:239)
DOYgroup <- factor(rep(1:3, times=c(3,3,3)))
DOYgroup <- factor(rep(c('171:173', '211:213', '237:239'), times=c(3,3,3)))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(bdat2014full$Time, na.rm=TRUE),
                                      max(bdat2014full$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
m3pred <- predict(m3, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, m3pred, DOYgroups)
names(predicted)[which(names(predicted)=='m3pred')] <- 'pH'
ggplot(predicted, aes(x = Time, y = pH, group= DOY, col=DOYgroups)) +
  theme_bw() +
  geom_line() +
  scale_colour_discrete(name="Day of Year") +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab('pH')

## for co2
N <- 200
simDOY <- c(168:170, 190:192, 228:230)
DOYgroup <- factor(rep(1:3, times=c(3,3,3)))
DOYgroup <- factor(rep(c('Jun 17-19', 'Jul 9-11', 'Aug 16-18'), times=c(3,3,3)))  
reptimes <- length(simDOY)
preddf <- data.frame(`Time` = rep(seq(min(bdat2014full$Time, na.rm=TRUE),
                                      max(bdat2014full$Time, na.rm=TRUE),
                                      length = N), times=reptimes),
                     `DOY` = rep(simDOY, each = N))
mco2pred <- predict(mco2, newdata = preddf, type = "link")

DOYgroups <- rep(DOYgroup, each = N)

predicted <- cbind(preddf, mco2pred, DOYgroups)
names(predicted)[which(names(predicted)=='mco2pred')] <- 'pCO2'

labdatco2 <- data.frame(x = 3, y = 380, label = "Mean (atm)")

co2modplot <- ggplot(predicted, aes(x = Time, y = pCO2, group= DOY, col=DOYgroups)) +
  theme_bw(base_size = 10) +
  scale_colour_brewer(name = "Date", type = 'qual', palette = 'Dark2', direction=1) +
  geom_line() +
  theme(legend.position="top") +
  xlab('Time of Day') + ylab(expression(paste(italic(p)*"CO"[2]~"(ppm)"))) +
  geom_abline(intercept = 398, slope = 0, linetype='dotted') +
  geom_text(data = labdatco2, aes(label = label, x = x, y = y), 
            show.legend = FALSE, inherit.aes = FALSE,  size = 3)
ggsave('../docs/private/bp-diel-co2gam.png', width=20, height=15, units='cm')

## arrange a few plots
asloplots <- plot_grid(co2modplot, dielplot, ncol=2,label_size = 10)
ggsave('../docs/private/dielplots.png', asloplots, width=20, height=10, units='cm')

## gam on other params
bdat2014fullf <- transform(bdat2014fullf, daylength = sunset-sunrise)
bdat2014fullf$daylength <- as.numeric(bdat2014fullf$daylength)
bpmod <- gam(ph1 ~
               s(chl) +
               s(airtemp) +
               s(daylength) +
               s(ODOrel1) +
               s(windsp),
             data = bdat2014fullf,
             select = TRUE, method = "REML", family = gaussian(),
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE,
                                   newton = list(maxHalf = 60)))
## FIXME: tested scat but though histogram looked better nothing else did...??
saveRDS(bpmod, '../data/private/bpmod.rds')

bpco2phmod <- gam(co2corr ~ s(ph1), method = "REML", family = gaussian,
                  na.action = na.exclude, data=bdat2014fullf,
                  control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
saveRDS(bpco2phmod, '../data/private/bpco2phmod.rds')
