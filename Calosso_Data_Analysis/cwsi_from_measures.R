rm(list=ls())

library(raster)
library(rgdal)
library(dplyr)
library(RColorBrewer)
source('CWSI_fun_noTd.R')


set.seed(1234)


# import the data
tab <- read.table('data_calosso/C.txt', sep = ';', header=T)
n <- nrow(tab)

plot(tab$CANOPY_TEMP, ylab = 'T_canopy_C', type = 'b')
grid()
#plot(tab$CANOPY_TEMP,type = 'b')


######################################
# outlier correction A
# tab$CANOPY_TEMP[2] <- 22
# tab$CANOPY_TEMP[31] <- 26
# tab$CANOPY_TEMP[37] <- 27.5
# tab$CANOPY_TEMP[10] <- 23
# tab$WIND_SPEED <- rep(mean(tab$WIND_SPEED), 42)


######################################
# outlier correction B
# tab$CANOPY_TEMP[42] <- mean(tab$CANOPY_TEMP, na.rm=T)
# tab$WIND_SPEED <- rep(mean(tab$WIND_SPEED, na.rm=T), 48)
# tab$UV_RAD[35] <- mean(tab$UV_RAD, na.rm=T)


######################################
# outlier correction C
tab$CANOPY_TEMP[1] <- 22
tab$CANOPY_TEMP[5] <- 22.5
tab$CANOPY_TEMP[13] <- 23
tab$CANOPY_TEMP[23] <- 24
tab$WIND_SPEED <- rep(mean(tab$WIND_SPEED, na.rm=T), 35)


###########################################################################################
# convert coordinates from Long-Lat to UTM
long_lat <- tab[,2:3]


names(long_lat) <- c("X","Y")

# Convert it to a sp object
coordinates(long_lat) <- ~ X + Y # longitude first

# Add a coordinate reference system
proj4string(long_lat) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Project using spTransform
Utm <- spTransform(long_lat, CRS("+proj=utm +zone=32 ellps=WGS84"))


xy <- Utm@coords




tab <- cbind(tab[,4:11], xy)

tab['x_grid'] <- NA
tab['y_grid'] <- NA


########################################################################################
# Define the grid

# set grid limits of the interested area
my_xmin <- 438765
my_xmax <- 439005
my_ymin <- 4954425
my_ymax <- 4954845
resolution <- 10
field_grid <- raster(extent(my_xmin, my_xmax, my_ymin, my_ymax), res=resolution, crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))

field_grid_tab <- rasterToPoints(field_grid)



########################################################################################
# find the cell in the field grid where each measurment is
for(i in 1:n){
  meas_x = tab$X[i]
  meas_y = tab$Y[i]
  
  min_dist_x <- min(abs(field_grid_tab[,1] - meas_x))
  min_dist_y <- min(abs(field_grid_tab[,2] - meas_y))
  
  meas_idx <- which(abs(field_grid_tab[,1] - meas_x) == min_dist_x & abs(field_grid_tab[,2] - meas_y) == min_dist_y)
  tab$x_grid[i] <- field_grid_tab[meas_idx,1]
  tab$y_grid[i] <- field_grid_tab[meas_idx,2]
}


# Original Measurment Points
plot(tab$X, tab$Y, xlim=c(my_xmin, my_xmax), ylim=c(my_ymin, my_ymax), main='Original Measurment Points')
abline(h=c(my_ymin, my_ymax), v=c(my_xmin, my_xmax), col='red')
abline(h=seq(my_ymin, my_ymax, resolution), v=seq(my_xmin, my_xmax, resolution), col='red')


# Adjusted Measurment Points
plot(tab$x_grid, tab$y_grid, xlim=c(my_xmin, my_xmax), ylim=c(my_ymin, my_ymax), main='Adjusted Measurment Points')
abline(h=c(my_ymin, my_ymax), v=c(my_xmin, my_xmax), col='red')
abline(h=seq(my_ymin, my_ymax, resolution), v=seq(my_xmin, my_xmax, resolution), col='red')



########################################################################################
# compute CWSI
# Ta_min <- 22 #A
# Ta_max <- 25 #A
# Ta_min <- 27.3  #B
# Ta_max <- 28.5  #B
Ta_min <- 25.0  #C
Ta_max <- 25.8  #C
#Ta_min <- 28  #D
#Ta_max <- 28  #D
Ta_corretto <- seq(Ta_min, Ta_max, length=n)
plot(Ta_corretto, type = 'b')

CWSI_1_5m <- rep(0,n)
for(i in 1:n){
  #CWSI_1_5m[i] <- CWSI_fun_noTd(Ta=(tab$AIR_TEMP[i]-3), Tc=tab$CANOPY_TEMP[i], UV=tab$UV_RAD[i], IR=tab$IR_RAD[i], p=tab$PRESSURE[i], Ur=tab$HUM[i], u2=tab$WIND_SPEED[i], h=1.5)
  CWSI_1_5m[i] <- CWSI_fun_noTd(Ta=Ta_corretto[i], Tc=tab$CANOPY_TEMP[i], UV=tab$UV_RAD[i],  p=tab$PRESSURE[i], Ur=tab$HUM[i], u2=tab$WIND_SPEED[i], h=1.5)
}


########################################################################################
# CWSI grid  (perform the mean in cells with more than one value)
CWSI_grid <- cbind(tab[,11:12], CWSI_1_5m)

CWSI_grid_no_dupl <- CWSI_grid %>% group_by(across(all_of(c('x_grid', 'y_grid')))) %>% summarise(avg=mean(CWSI_1_5m))



#######################################################################################
# create a raster
raster_cwsi <- rasterFromXYZ(CWSI_grid_no_dupl, crs='+proj=utm +zone=32 ellps=WGS84')

plot(raster_cwsi, col=rev(brewer.pal(10, 'RdYlGn')), main='CWSI 1.5m', alpha=0.75)

## griglie dai dati rilevati 
# writeRaster(raster_cwsi,'CWSI_from_meas_A.tif',format="GTiff", overwrite=T)
# write.table(CWSI_grid, 'CWSI_from_meas_A.csv', row.names = F)

## griglie con i dati sistemati - outlier corrections
writeRaster(raster_cwsi,'CWSI_corrected_C.tif',format="GTiff", overwrite=T)
write.table(CWSI_grid, 'CWSI_corrected_C.csv', row.names = F)


#######################################################################################
#Ground humidity raster
GH_grid <- cbind(tab[,11:12], tab$GROUND_HUM)

#GH_grid_no_dupl <- GH_grid %>% group_by(across(all_of(c('x_grid', 'y_grid')))) %>% summarise(avg=mean(tab$GROUND_HUM))

raster_GH <- rasterFromXYZ(GH_grid, crs='+proj=utm +zone=32 ellps=WGS84')

plot(raster_GH, col=rev(brewer.pal(10, 'RdYlGn')), main='Ground humidity', alpha=0.75)

cor(tab$GROUND_HUM,CWSI_1_5m, use='pairwise.complete.obs')

#writeRaster(raster_GH,'GH_from_meas_A.tif',format="GTiff", overwrite=T)



#cwsi_points
CWSI_points <- cbind(tab[,9:10], CWSI_1_5m)
write.csv(CWSI_points, 'CWSI_points_A.csv')

        



