.libPaths("~/R")
setwd("~/Documents/R/Clark")

#libraries

require(raster)
require(dplyr)
require(spData)
require(ggplot2)
require(maps)
require(gstat)

#read in data 

pr <- read.csv("pr_yearly.csv", stringsAsFactors = F)
eto <- read.csv("eto_yearly.csv", stringsAsFactors = F)
tmin <- read.csv("tmin_yearly.csv", stringsAsFactors = F) #make sure to rename columns from this file
col <- which(names(tmin) == "eto")
colnames(tmin)[col] <- "tmin"
tmin$tmin <- tmin$tmin - 273.5 #convert to celcius

#create blank raster 

proj <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

points <- SpatialPointsDataFrame(coords= tmin[, c("LON", "LAT")], 
                                 proj4string = CRS(proj), 
                                 data = tmin)

blankRaster <- raster(ncol = as.integer(2 * round(max(tmin$LON) - min(tmin$LON), 0)),
                                        nrow = as.integer(2 * round(max(tmin$LAT) - min(tmin$LAT)), 0), 
                                        xmn = min(tmin$LON), xmx = max(tmin$LON), ymn = min(tmin$LAT), 
                                        ymx = max(tmin$LAT), crs = CRS(proj))

r1 <- rasterize(points, blankRaster, field = points$tmin, fun = mean, background = NA, na.rm = T)
plot(r1) #plots mean of tmin from 1995-2018


#create ts model for temperature change 
agg.tmin <- aggregate(tmin ~ year, data = tmin, FUN = mean, na.rm = T)
tmin.mod <- lm(tmin ~ year, data = agg.tmin)


#visualize 

plot(agg.tmin$year, agg.tmin$tmin, type = "l")
abline(lm(agg.tmin$tmin ~ agg.tmin$year), col = "red")

temporalG <- tmin.mod$coefficients[2]

#calculating spatial gradient 


minYear <-min(unique(tmin$year))
maxYear <-max(unique(tmin$year))

minY <- tmin[tmin$year == minYear, ]
maxY <- tmin[tmin$year == maxYear, ]

changeT <- maxY$tmin - minY$tmin

changeT <- data.frame(cbind(minY[, -c(which(names(minY) == "tmin"))], changeT))
changeT <- changeT[, -1]

#create raster for change in min temp from 1995 to 2018 (in C)
proj <- "+proj=aea +lat_1=20 +lat_2=6 +lat_0=40 +lat_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83"

points <- SpatialPointsDataFrame(coords= changeT[, c("LON", "LAT")], 
                                 proj4string = CRS(proj), 
                                 data = changeT)

blankRaster <- raster(ncol = as.integer(2 * round(max(changeT$LON) - min(changeT$LON), 0)),
                                        nrow = as.integer(2 * round(max(changeT$LAT) - min(changeT$LAT)), 0), 
                                        xmn = min(changeT$LON), xmx = max(changeT$LON), ymn = min(changeT$LAT), 
                                        ymx = max(changeT$LAT), crs = CRS(proj))

r1 <- rasterize(points, blankRaster, field = points$changeT, fun = mean, background = NA, na.rm = T)

plot(r1) #plots change in minimum temperature in C

#spatial gradient functions 


get_NS_diffs <- function(avg_raster) {
  f <- matrix(rep(1, 9), nrow=3, ncol=3)
  NS_diff <- raster::focal(avg_raster, w = f, nrow=3, ncol=3, pad = TRUE,
                           fun = function(x, ...) {
                             be <- x[2] - x[5]
                           }
  )
  NS_adj_diff <- NS_diff / 111.325
  return(NS_adj_diff)
}

NS_gradient <- get_NS_diffs(r1)


get_WE_diffs <- function(avg_raster) {
  f <- matrix(rep(1, 9), nrow=3, ncol=3)
  lat <- sp::coordinates(avg_raster)[, 2]
  WE_diff <- raster::focal(avg_raster, w = f, nrow=3, ncol=3, pad = TRUE,
                           fun = function(x, ...) {
                             ba <- x[6] - x[5]
                           }
  )
  WE_adj_diff <- WE_diff / (111.325 * cos (lat * pi / 180))
  return(WE_adj_diff)
}

WE_gradient <- get_WE_diffs(r1)


get_spatial_gradient <- function(NS_gradient, WE_gradient) {
  f <- matrix(rep(1, 9), nrow=3, ncol=3)
  x_gradient <- raster::focal(WE_gradient, w = f, nrow=3, ncol=3, pad = TRUE,
                              fun = function(x, ...) {
                                #browser()
                                wt <- c(1, 2, 1, 2, 0, 2, 1, 2, 1)
                                weighted.mean(x, wt, na.rm = TRUE)
                              })
  y_gradient <- raster::focal(NS_gradient, w = f, nrow = 3, ncol = 3, pad = TRUE,
                              fun = function(x, ...) {
                                wt <- c(1, 2, 1, 2, 0, 2, 1, 2, 1)
                                weighted.mean(x, wt, na.rm = TRUE)
                              })
  
  # Get magnitude of resultant vector
  magnitude <- sqrt(x_gradient^2 + y_gradient^2)
  
  # Get angle of resultant vector
  angle <- raster::atan2(x_gradient, y_gradient) * 180 / pi
  
  # Create correction raster to produce positive angles
  neg_angles <- (angle < 0) * 360
  
  angle <-
    raster::overlay(x = angle, y = neg_angles, fun = function(x, y) {x + y})
  
  spatial_grad_brick <- raster::brick(magnitude, angle)
  names(spatial_grad_brick) <- c('spatial_gradient', 'angle')
  return(spatial_grad_brick)
}

spatialG <- get_spatial_gradient(NS_gradient, WE_gradient)
spatialG1 <-spatialG*100

#temporal gradient function 




get_sst_linear_change <- function(hadsst_raster = b, years = 2000:2011) {
  # annual_rasters <- get_annual_ssts(hadsst_raster, years)
  
  time_ <- I(years - mean(years))
  fun <- function(x) {
    # if (sum(is.na(x)) < length(x)) {
      slope <- 10 * lm(x ~ time_, na.action = na.exclude)$coefficients[2]
      return(slope)
    # }
    # return(NA)
  }
  lin_change <- raster::calc(b, fun)
  return(lin_change)
}

tminAvgRaster <- list()
year <- unique(tmin$year)
proj <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"


for(i in 1:length(year)){
  
  tmp <- tmin[tmin$year == year[i], ]
  points <- SpatialPointsDataFrame(coords= tmp[, c("LON", "LAT")], 
                                   proj4string = CRS(proj), 
                                   data = tmp)
  
  blankRaster <- raster(ncol = as.integer(2 * round(max(tmp$LON) - min(tmp$LON), 0)),
                        nrow = as.integer(2 * round(max(tmp$LAT) - min(tmp$LAT)), 0), 
                        xmn = min(tmp$LON), xmx = max(tmp$LON), ymn = min(tmp$LAT), 
                        ymx = max(tmp$LAT), crs = CRS(proj))
  
  tminAvgRaster[[i]] <- rasterize(points, blankRaster, field = points$tmin, fun = mean, background = NA, na.rm = T)
  
  
}

tminAvg <- brick(tminAvgRaster)
tminAvg1 <- mean(tminAvgRaster)

temporalG <- get_sst_linear_change(tminAvg[1], years = 1995:2018)
years <- 1995:2018
time <- I(years - mean(years))

slopeFun <- function(tmin, time){
  slope <- 10 * lm(tmin ~ time, na.action = na.exclude)$coefficients[2]
return(slope)
}

t <- raster::calc(tminAvg, trial)
geeIDlength <- (unique(tmin$geeID))
slopeAll <- matrix(NA, nrow = length(geeIDlength), ncol = 2)


for(i in 1:length(geeIDlength)){
  
  tmp <- tmin[tmin$geeID == geeIDlength[i], ]
  
  tmp1 <- slopeFun(tmp$tmin, time =years)
  
  slopeAll[i, 1] <- tmp1
  slopeAll[i, 2] <- geeIDlength[i]
  
  if(i %% 1000 == 0){print(paste0((304140 - i), " iterations left!"))}
  
}

colnames(slopeAll) <- c("TemporalG", "geeID")
slopeAll <- as.data.frame(slopeAll)
slopeAll1 <- left_join(slopeAll, tmin, by = "geeID")


points <- SpatialPointsDataFrame(coords= slopeAll1[, c("LON", "LAT")], 
                                 proj4string = CRS(proj), 
                                 data = slopeAll1)

blankRaster <- raster(ncol = as.integer(2 * round(max(slopeAll1$LON) - min(slopeAll1$LON), 0)),
                      nrow = as.integer(2 * round(max(slopeAll1$LAT) - min(slopeAll1$LAT)), 0), 
                      xmn = min(slopeAll1$LON), xmx = max(slopeAll1$LON), ymn = min(slopeAll1$LAT), 
                      ymx = max(slopeAll1$LAT), crs = CRS(proj))

temporalG <- rasterize(points, blankRaster, field = (points$TemporalG), fun = mean, background = NA, na.rm = T)

plot(temporalG)
plot(spatialG[[1]])
linear_change / raster::subset(spatial_gradient, 'spatial_gradient')
velocity <- temporalG/raster::subset(spatialG1, 'spatial_gradient')
plot(velocity)

