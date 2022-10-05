rm(list=ls())

library(raster)
library(rgdal)
library(dplyr)
source('CWSI_fun.R')
source('functions.R')


seed(1234)


# import the data
tab <- read.table('tab.txt')
n <- nrow(tab)


###########################################################################################
# convert coordinates from Long-Lat to UTM
long_lat <- tab[,9:10]


names(long_lat) <- c("X","Y")

# Convert it to a sp object
coordinates(long_lat) <- ~ X + Y # longitude first

# Add a coordinate reference system
proj4string(long_lat) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# Project using spTransform
Utm <- spTransform(long_lat, CRS("+proj=utm +zone=32 ellps=WGS84"))


xy <- Utm@coords


# Add radom in x y coordinates
xy[,1] <- xy[,1] + rnorm(n, mean=0, sd=6)
xy[,2] <- xy[,2] + rnorm(n, mean=0, sd=6)


tab <- cbind(tab[,1:8], xy)

tab['x_grid'] <- NA
tab['y_grid'] <- NA


########################################################################################
# Define the grid

# set grid limits of the intersted area
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

t <- 5

plot(tab$X, tab$Y)
abline(h=c(my_ymin, my_ymax), v=c(my_xmin, my_xmax), col='red')
points(tab[1:t, 9:10], col=rainbow(t), pch=19)
abline(h=seq(my_ymin, my_ymax, resolution), v=seq(my_xmin, my_xmax, resolution), col='red')


plot(tab$x_grid, tab$y_grid)
abline(h=c(my_ymin, my_ymax), v=c(my_xmin, my_xmax), col='red')
points(tab[1:t, 11:12], col=rainbow(t), pch=19)
abline(h=seq(my_ymin, my_ymax, resolution), v=seq(my_xmin, my_xmax, resolution), col='red')


########################################################################################
# compute CWSI
CWSI_1_5m <- rep(0,n)
for(i in 1:n){
  CWSI_1_5m[i] <- CWSI_fun(tab$Ta[i], tab$Tc[i] + 273.15, tab$Rn[i]/3600, tab$p[i]/100, tab$Td[i], tab$ux[i], tab$uy[i], h=1.5)
}


########################################################################################
# CWSI grid  (perform the mean in cells with more than one value)
CWSI_grid <- cbind(tab[,11:12], CWSI_1_5m)

CWSI_grid_no_dupl <- CWSI_grid %>% group_by(across(all_of(c('x_grid', 'y_grid')))) %>% summarise(avg=mean(CWSI_1_5m))



#######################################################################################
# create a raster
raster_cwsi <- rasterFromXYZ(CWSI_grid_no_dupl, crs='+proj=utm +zone=32 ellps=WGS84')

plot(raster_cwsi)

writeRaster(raster_cwsi,'CWSI_from_meas.tif',format="GTiff", overwrite=T)
write.table(CWSI_grid, 'CWSI_from_meas.csv', row.names = F)
