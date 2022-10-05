rm(list=ls())
library(dplyr)
library(lubridate)
library(tidyr)
library(robustbase)
library(psych)
library(MASS)
library(ellipse)
library(here)
library(DescTools)
library(knitr)
library(RobStatTM)


tab <- read.csv('CWSI15m.csv')
tab$date <- as.Date(tab$date, "%Y-%m-%d")


tab$month <- format(tab$date,'%m')


# keep data only between March and September
tab <- tab %>% filter(month %in% c("03", "04", "05", "06", "07", "08", "09"))
nn <- dim(tab)[1]


# Assign season Spring from March to June and Summer from July to September
tab$season <- rep(0, nn)
tab$season[tab$month %in% c('03', '04','05', '06')] <- 'spring'
tab$season[tab$month %in% c('07','08','09')] <- 'summer'


# compute CWSI mean by year and season
CWSI_means <- tab %>% group_by(year, season) %>% summarise(CSWI_mean = mean(CWSI_1_5m, na.rm=TRUE))

# compute CWSI standard deviation by year
CWSI_sd <- tab %>% group_by(year) %>% summarize(CWSI_sd = sd(CWSI_1_5m, na.rm=TRUE))


# create the table
tab_means <- CWSI_means %>% spread(season, CSWI_mean)


# add productivity
prods <- c(270, 425, 254, 255, 552, 547, 466, 592, 327, 595, 610, 437, 580, 591, 432, 468, 595, 370, 598, 461, 513, 604)
tab_prod <- data.frame(year=tab_means$year, prod=prods)

tab_data <- inner_join(tab_means, tab_prod, by='year')


# remove old years
tab_data <- tab_data %>% filter(year >= 2010 )



# Spring
plot(tab_data$spring, tab_data$prod, pch=19)

alph <- 0.95     # percentage of data to use
mod_spring <- ltsReg(prod ~ spring, alpha=alph, mcd=TRUE, data=tab_data) 
summary(mod_spring)

plot(tab_data$spring, tab_data$prod, pch=19,
     xlab='CWSI Spring', ylab='Production [100kg]')
abline(mod_spring, col='red', lwd=2)



# Summer
plot(tab_data$summer, tab_data$prod, pch=19)

alph <- 0.95    # percentage of data to use
mod_summer <- ltsReg(prod ~ summer, alpha=alph, mcd=TRUE, data=tab_data) 
summary(mod_summer)

plot(tab_data$summer, tab_data$prod, pch=19,
     xlab='CWSI Summer', ylab='Production [100kg]', main='Regression Model')
abline(mod_summer, col='red', lwd=2)



# Spring and Summer
alph <- 0.95     # percentage of data to use
Z <- tab_data[,2:4]
mod_both <- ltsReg(prod ~ ., alpha=alph, mcd=TRUE, data=Z) 
summary(mod_both)


# diagnostic
plot(mod_both)



#############################################
# Plot IC
xx = seq(0.45, 0.75, length=1000)
X = data.frame(summer=xx)

Z_rm = Z[-5,]
mod_lm = lm(prod ~ summer,  data=Z_rm)

IC_pred = predict(mod_lm, X, interval='confidence')

plot(tab_data$summer, tab_data$prod, pch=19, type='n', xlab='CWSI Summer', ylab='Production [100kg]', main='Regression Model')
polygon(c(xx, rev(xx)), c(IC_pred[,3], rev(IC_pred[,2])), col = "gray")
lines(xx, IC_pred[,1], col='red', lwd=2)
#abline(mod_lm, col='red', lwd=2)
points(tab_data$summer, tab_data$prod, pch=19)



# Spring and Summer Interaction
alph <- 0.95     # percentage of data to use
Z <- tab_data[,2:4]
mod_interact <- ltsReg(prod ~ spring:summer, alpha=alph, mcd=TRUE, data=Z) 
summary(mod_interact)





