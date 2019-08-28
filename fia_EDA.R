.libPaths("/Users/Shubhi/Documents/R")
setwd("/Users/Shubhi/Documents/R/Clark")
require(sf)
require(raster)
require(dplyr)
require(spData)
# require(spDataLarge)
require(ggplot2)
require(maps)

t <- read.csv("WY_AllPlots.csv") 

allFIA <- read.csv("FIA_150419.csv")
allFIA$BA <- (pi*(allFIA$DIA/2)^2)/allFIA$TPA_UNADJ
WY <- allFIA[allFIA$STATE == "WY",]

allFIA$g.species <- paste(allFIA$genus, "-", allFIA$species)

acer.rubrum <- allFIA[allFIA$g.species == "Acer - rubrum",]

usa <- map_data("usa")
states <- map_data("state")
ggplot() + geom_polygon(data = usa, aes(x=long, y = lat, group = group)) + 
  coord_fixed(1.3)
#acer.rubrum ---- 
ggplot(acer.rubrum, aes(x = acer.rubrum$LON, y = acer.rubrum$LAT))+ geom_point(color = "green", alpha = 1/10) + theme_bw() + ggtitle("Acer Rubrum (red maple)") +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3) 
#pinus.taeda ----
pinus.taeda <- allFIA[allFIA$g.species == "Pinus - taeda", ]

ggplot(pinus.taeda, aes(x = pinus.taeda$LON, y = pinus.taeda$LAT)) + 
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3) + geom_point(color = "green", alpha = 1/10) + theme_bw() + ggtitle("Loblolly")


# Liquidambar styraciflua -----
liqui.stry <- allFIA[allFIA$g.species == "Liquidambar - styraciflua", ]
ggplot(liqui.stry, aes(x = liqui.stry$LON, y = liqui.stry$LAT)) + geom_point(color = "green", alpha = 1/10) + theme_bw() + ggtitle("Sweetgum") +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3) 

# Pseudotsuga menziesii-----

pseudo.menzi <- allFIA[allFIA$g.species == "Pseudotsuga - menziesii", ]

ggplot(pseudo.menzi, aes(x = pseudo.menzi$LON, y = pseudo.menzi$LAT)) + geom_point(color = "green", alpha = 1/10) + theme_bw() + ggtitle("Douglas Fir") +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3) 

# Populus tremuloides-----

pop.tremul <- allFIA[allFIA$g.species == "Populus - tremuloides",]

ggplot(pop.tremul, aes(x = pop.tremul$LON, y = pop.tremul$LAT)) + geom_point(color = "green", alpha = 1/10) + theme_bw() + ggtitle("Quaking Aspen") +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3) 

# Abies balsamea  -----

abies.balsamea <- allFIA[allFIA$g.species == "Abies - balsamea", ]

ggplot(abies.balsamea, aes(x = abies.balsamea$LON, y = abies.balsamea$LAT)) + geom_point(color = "green", alpha = 1/10) + theme_bw() + ggtitle("Balsam Fir") +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3) 

#Ponderosa pine ----- 

pinus.ponderosa <- allFIA[allFIA$g.species == "Pinus - ponderosa", ]

ggplot(pinus.ponderosa, aes(x = pinus.ponderosa$LON, y = pinus.ponderosa$LAT)) + geom_point(color = "green", alpha = 1/10) + theme_bw() + ggtitle("Pinus ponderosa") +
  geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
  coord_fixed(1.3)

t <- pinus.ponderosa[pinus.ponderosa$STATE == "WY",]

#latest FIA EDA ------

fia <- read.csv("FIA_061219.csv")
fiaDBF <- read.dbf("~/Documents/R/Clark/FIA/shape/fia.dbf")

uniqueDBF <- unique(fiaDBF$cood)
uniqueFIA <- unique(fia$cood)





#remove trees if STATUSCD = 2 or STATUSCD = 3
#potentially remove plots if STATUSCD = 2 or STATUSCD = 3

fia1 <- fia[-c(which(fia$STATUSCD == 2 | fia$STATUSCD == 3)), ]
fia1 <- fia1[-which(fia1$TPA_UNADJ > 100), ]
fia1 <- fia1[-c(which(fia1$PCTFOREST < 0.25)), ]

TPAna <- which(is.na(fia1$TPA_UNADJ))
PCTna <- which(is.na(fia1$PCTFOREST))

fia1 <- fia1[-c(TPAna, PCTna), ]

fia1$BA <- ((pi*(fia1$DIA/2)^2)*fia1$TPA_UNADJ)/fia1$PCTFOREST


# ba.agg <- aggregate(fia1$BAtmp, by = list(fia1$IDplot, fia1$INVYR, fia1$cood), sum, na.rm = T)

ba.agg <- aggregate(BAtmp ~ IDplot + INVYR + cood + LAT + LON, data = fia1, FUN= sum, na.rm= T)
#fia plots <- remove duplicated lat lon- plotID 

ba.agg$ba.agg.metric <-  (ba.agg$BAtmp*0.00064516)*2.47

# ba.agg$ba.agg.wo <- ba.agg.metric[which(ba.agg.metric <1000)]


# fia.ba <- fia[-c(which(fia$STATUSCD == 2 | fia$STATUSCD == 3)), ]

# BA <- (((pi*(dia/2)^2)* fia$TPA_UNADJ)/fia$PCTFOREST)/2.47 #10-100s 

# BA.df <- data.frame(cbind(BA = fia$BA, LAT = fia$LAT, LON = fia$LON, STATE = fia$STATE)) 

# ggplot(ba.agg, aes(x = ba.agg$LON, y = ba.agg$LAT)) + geom_polygon(data = usa, aes(x=long, y = lat, group = group), fill = NA, color = "black") + geom_point(aes(fill = ba.agg$ba.agg.metric)) + theme_bw()+ 
  # coord_fixed(1.3)



####### fia file for upload 
xvar <- c(which(colnames(fia) =="ASPECT"),which(colnames(fia) =="cood"),
          which(colnames(fia) == "ELEV"),which(colnames(fia) == "STDAGE"),
          which(colnames(fia) == "SLOPE"))
fia.xvar <- fia[, (xvar)]

write.csv(fia.xvar, "fia.xvar.csv")

dropCols <- c(which(names(ba.agg) == "IDplot"),which(names(ba.agg) == "INVYR"), which(names(ba.agg) == "LAT"), 
                  which(names(ba.agg) == "LON"))

ba.agg.f <- ba.agg[, -dropCols]

write.csv(ba.agg.f, "fiaBA.csv")

###### FIA y data for upload 

# ba.agg$ba.agg.metric <-  (ba.agg$BAtmp*0.00064516)*2.47

fia1$BA <- (fia1$BA * 0.00064516)*2.47

# ba.agg <- aggregate(BAtmp ~ IDplot + INVYR + cood + LAT + LON, data = fia1, FUN= sum, na.rm= T)

BA <- aggregate(BA ~ IDplot + INVYR + cood + SPCD, data = fia1, FUN = sum, na.rm= T)

write.csv(BA, "BA.csv")

#create raster 

proj <- "+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"

points <- SpatialPointsDataFrame(coords = ba.agg[, c("LON", "LAT")], 
                                 proj4string = CRS(proj), 
                                 data = ba.agg)

blankRaster <- raster(ncol = as.integer(2*round((max(ba.agg$LON) - min(ba.agg$LON)), 0)), nrow = as.integer(2*round((max(ba.agg$LAT) - min(ba.agg$LAT)), 0)),
                      xmn = min(ba.agg$LON), xmx = max(ba.agg$LON), ymn = min(ba.agg$LAT), ymx = max(ba.agg$LAT), 
                      crs = CRS(proj))

r1 <- rasterize(points, blankRaster, field = points$ba.agg.metric, fun = mean, background = NA, na.rm = T)
