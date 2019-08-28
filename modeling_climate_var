# script to cluster design matrix data, 
# cluster y matrix data and run GJAM on entire FIA dataset 

#libraries -----

require(dplyr)
require(gjam)
require(foreign)
source('clarkFunctions.R')


#functions ----

#organize climate variable (returns average annual)
FIAClimVar <- function(csvPath, dbfPath, climVars, minYear, maxYear){
  
  csvPath <- paste0(getwd())
  # dbfPath <- paste0(getwd(), "/FIA/shape/fia.dbf")
  # minYear = minYear
  # maxYear = maxYear
  # climVars = "pr"
  
  # fia <- read.dbf(file = dbfPath)
  clim <- read.csv(paste0(csvPath, "/", climVars, "_monthly_G.csv"))
  
  clim$year <- substr(clim$yr_m, 1, 4)
  
  clim <- clim[, -(which(names(clim) == "yr_m"))]
  clim <- filter(clim, clim$year >= minYear)
  
  
  col <- which(colnames(clim) == climVars)
  byCol1 <- which(colnames(clim) == "geeID" )
  byCol2 <- which(colnames(clim) == "year")
  clim.agg <- aggregate(clim[,col], by = list(clim[, byCol1], clim[, byCol2]), FUN= mean, na.rm= TRUE)
  colnames(clim.agg) <- c("geeID", "year", climVar)
  clim.all <- merge(fia, clim.agg, by = "geeID", all.y = F, sort = T)
  return(clim.all)
}


#trim dataset
subsetCols <- function(varList, dataset){
  cols <- vector(length = length(varList))
  for(i in 1:length(varList)){
    cols[i] <- which(names(dataset) == varList[i])
  }
  return(cols)
}


#creating yearly csv

fia <- read.dbf("fia.dbf")
fia$cood <- paste0(fia$LAT, "_", fia$LON)

pr <- FIAClimVar(csvPath =clim.path, climVars = "pr", minYear = 1995, maxYear = 2018)
# write.csv(pr, "pr_yearly.csv")
eto <- FIAClimVar(csvPath = csvPath, climVars = "eto", minYear = 1995, maxYear = 2018)
# write.csv(eto, "eto_yearly.csv")
pr_eto <- pr$pr - eto$eto
rawPrEto <- cbind(pr, pr_eto, eto$eto)

tmin <- FIAClimVar(csvPath = csvPath, climVars = "tmin", minYear = 1995, maxYear = 2018)
# write.csv(tmin, "tmin_yearly.csv")

#once files are created, read them in as,
#tmin <- read.csv("tmin_yearly.csv")
#colnames(tmin)[which(names(tmin) == "eto")] <- "tmin"
#eto <- read.csv("eto_yearly.csv")
#pr <- read.csv("pr_yearly.csv")

prEto <- pr$pr - eto$eto
prEto.df <- data.frame(cbind(pr, prEto))
prEto <- prEto.df[, c(subsetCols(c("geeID", "year", "LON", "LAT", "prEto"), prEto.df))]
prEto$cood <- paste0(prEto$LAT, "_", prEto$LON)
prEto$year <- as.numeric(prEto$year)

tmin$cood <- paste0(tmin$LAT, "_", tmin$LON)
# incorrectCol <- which(names(tmin) == "eto")
# colnames(tmin)[incorrectCol] <- "tmin"



climVar <- left_join(prEto, tmin, by = c("geeID", "year"))
climVar <- climVar[, subsetCols(c("geeID", "year", "LON.x", "LAT.x", 
                                  "prEto", "cood.x", "tmin"), climVar)]

FIA<- read.csv("FIA_061219.csv")

FIAx <- FIA[, subsetCols(c("ASPECT", "INVYR", "ELEV", "STDAGE",
                              "SLOPE", "LON", "LAT","IDplot"), FIA)]

# rm(FIA)
FIAx <- rbind(FIAx[!(is.na(FIAx$ELEV)), ] ,FIAx[(is.na(FIAx$ELEV)), ])

FIAx <- FIAx[!(duplicated(FIAx[, c("IDplot", "INVYR")])), ]

FIAx$cood <- paste0(FIAx$LAT, "_", FIAx$LON)

soil <- read.csv("SoilsFIA.csv")
soil <- soil[, subsetCols(c("LAT", "LON", "geeID", "clay", "sand", "soilDepth"), soil)]
soil$cood <- paste0(soil$LAT, "_", soil$LON)


climSoil <- left_join(tabx4, soil, by = c("geeID"))
# climSoil <- climSoil[, c(subsetCols(c("geeID", "year", "LON", "LAT", 
                                      # "prEto", "cood", "tmin", "clay", "sand", "soilDepth"), climSoil))]


allX <- merge(FIAx, climSoil, by.x = c("cood", "INVYR"), by.y = c("cood", "year"), all.y = F)

allX <- allX[, c(subsetCols(c("cood", "IDplot","INVYR", "ASPECT", "ELEV", "SLOPE", 
                             "LON.x", "LAT.x", "tmin ,"clay", "sand",
                              "soilDepth", "prEto"), allX))]
cols <- c(which(names(allX) == "LON.x"), which(names(allX) == "LAT.x"))
colnames(allX)[cols] <- c("LON", "LAT")

allX$u1 <- sin(allX$SLOPE)
allX$u2 <- sin(allX$SLOPE)*sin(allX$ASPECT)
allX$u3 <- sin(allX$SLOPE)*cos(allX$ASPECT)

dropCols <- c(which(names(allX) == "SLOPE"), which(names(allX) == "ASPECT"))

allX <- allX[, -dropCols]

# LOOP 
#cluster every year 

year <- unique(allX$INVYR)
year <- year[order(year)]
# dir.create("~/clusterFiles/")
# dir.create("~/clusterXvar/")

for(i in 1:length(year)){
  
  print(paste("Starting iteration", i, "..."))
  allX.sub <- allX[allX$INVYR == year[i], ]
  allX.sub <- na.omit(allX.sub)
  
  if(nrow(allX.sub) == 0){
    print("Error: No rows in allx.sub variable after na.omit()")
 
    next
    }
  
  cluster <- allX.sub[, c(subsetCols(c("u1", "ELEV","prEto", "tmin" ,"clay", "sand",
                                       "soilDepth"), allX.sub))]
  
  cluster1 <- scale(cluster)
  cluster1[is.na(cluster1)] <- 0
  
  print(paste("Starting clustering..."))
  
  xCluster <- kmeansEqualGroups(cluster1, 13, 15, progress = F)
  
  allX.sub$cluster <- xCluster$cluster
  
  #csv plotID, INVYR, cood and clusterID
  tmp <- allX.sub[, c(subsetCols(c("cood", "INVYR", "IDplot","cluster"), allX.sub))]
  
  
  write.csv(tmp, paste0("~/clusterFiles/clusterID_", year[i], ".csv"))
  
  print("Aggregating...")
  cluster$cluster <- xCluster$cluster
  t <- melt(cluster, id = "cluster")
  agg <- dcast(t, cluster ~ variable, mean, na.rm= T)
  
  print("Finalizing data frame...")
  
  join <- left_join(allX.sub, agg, by = "cluster")
  
  join <- join[, c(subsetCols(c("cood", "INVYR", "IDplot.y", "cluster", "ELEV.y", 
                                "clay.y", "sand.y", 
                                "soilDepth.y", "u1.y", "prEto.y", "tmin.y"), join))]
  
  colnames(join) <- c("cood","INVYR", "IDplot","cluster", "ELEV",
                      "clay", "sand", "soilDepth", "u1", "prEto", "tmin")
  
  write.csv(join, paste0("~/clusterXvar/clusterXvar_", year[i], ".csv"))
  
  print(paste0("Finished iteration ", i ))
  
}

input.path <- "~/clusterXvar/"
clusterFiles <- dir(input.path)
allCluster <- list()
yearW <- list()
for( i in 1:length(clusterFiles)){
  
  allCluster[[i]] <- read.csv(paste0(input.path, clusterFiles[i]), stringsAsFactors = F)

  rg <- regexpr("_", clusterFiles[i])
  yearW[i] <- substr(clusterFiles[i], rg+1, rg+4)
  }

cluster.final <- do.call(rbind, allCluster)
dropCols <- c(which(names(cluster.final)== "cood"), which(names(cluster.final)== "IDplot"),
              which(names(cluster.final)== "X"))
cluster.final1 <- cluster.final[, -dropCols]
d <- duplicated(cluster.final1)
cluster.final1 <- cluster.final1[which(d == FALSE), ]

yearW <- do.call(rbind, yearW)

#Loop for BA aggregation 

baRaw <- read.csv("BA.csv")
ba <- baRaw[, -1]

ba <- ba[ba$INVYR %in% year, ]

dir.create("~/clusterBA/")
input.path <- "~/clusterFiles/"
clusterFiles <- dir(input.path)
allBA <- list()

yearW <- as.numeric(yearW)

for(i in 1:length(yearW)){
  print(paste("Starting iteration", i, "..."))
  
  cid <- read.csv(paste0(input.path, clusterFiles[i]), stringsAsFactors = F)
  cid <- cid[, -1]
  
  if(nrow(cid) == 0){
    print(paste0("Clustering for year ", yearW[i], " not done"))
    next
  }
  
  
  BA.sub <- ba[(ba$INVYR == yearW[i]), ]
  BA.sub$cood <- as.character(BA.sub$cood)
  
  all <- left_join(BA.sub, cid, by = c("IDplot", "INVYR"))
  all <- na.omit(all)
  # col <- which(names(all) == "INVYR.x")
  # colnames(all)[col] <- "INVYR"
  
  print("Aggregating...")
  
  agg <- aggregate(BA ~ IDplot + SPCD + INVYR + cluster, data = all, FUN = mean, na.rm = T )
  
  print(paste0("No.of unique species in year ", yearW[i], " is ", length(unique(agg$SPCD))))
  write.csv(agg, paste0("~/clusterBA/BA_", yearW[i], ".csv"))
  
  print(paste("Finished with iteration", i))
  
  allBA[[i]] <- agg
}


#binding BA files 

# input.path <- "~/clusterBA/"
# clusterBA <- dir(input.path)
# allBA <- list()
# for( i in 1:length(clusterBA)){
#   
#   allBA[[i]] <- read.csv(paste0(input.path, clusterBA[i]), stringsAsFactors = F)
# }


cluster.ba <- bind_rows(allBA)

sp.info<- read.csv('REF_SPECIES.csv',header=T)
mm<- match(cluster.ba$SPCD,sp.info$SPCD)
genus<- sp.info$GENUS[mm]
species<- sp.info$SPECIES[mm]


symbol <- sp.info$SPECIES_SYMBOL[mm]

cluster.ba$genus<- genus
cluster.ba$species<- species
cluster.ba$SPsymbol <- symbol

ba.sub <- cluster.ba[, c(subsetCols(c("BA", "SPsymbol","INVYR", "cluster"),
                                    cluster.ba))]

melt.ba <- melt(ba.sub, id = c("SPsymbol", "INVYR",  "cluster"))
ba.cast <- reshape2::dcast(melt.ba, INVYR + cluster+ variable ~ SPsymbol, fun.aggregate = mean, na.rm= T)
ba.cast[is.na(ba.cast)] <- 0

c <- colSums(ba.cast[, 4:ncol(ba.cast)])
c <- c[order(c, decreasing =  T)]

cnames <- names(c[1:100])

cnames <- c("INVYR",  "cluster", cnames)
ba.final  <- ba.cast[, names(ba.cast) %in% cnames]

#merging dataframes

xy <- merge(ba.final, cluster.final1, by = c("cluster", "INVYR"), all.y = F)
# write.csv(xy, "xydata_fiaGJAM_run2.csv")

y <- xy[, 3:102]
x <- xy[, 103:ncol(xy)]

#sample 70% randomly, save 30% as OOB samples

set.seed(072)
sample <- sample(seq(1, nrow(x)), 0.7*nrow(x))

y.gjam <- y[sample, ]
x.gjam <- x[sample, ]


#OOB samples 

oob <- which(!(seq(1, nrow(x)) %in%sample))
x.oob <- x[oob,]
y.oob <- y[oob,]
#write.csv(oob, "OOBSerialNo_run2.csv")
# oob <- read.csv("OOBSerialNo.csv")
# oob <- as.numeric(oob$x)

#GJAM

# form <- as.formula("~ prEto + tmin + clay + sand+ soilDepth + u1+ u2 + u3 + ELEV + STDAGE") - fiaGJAM2
# form <- as.formula("~ prEto + tmin + clay + sand + soilDepth + u1 + ELEV + STDAGE") - fiaGJAM_run2
# form <- as.formula("~prEto + tmin + clay + sand + soilDepth + u1 + ELEV") 
# - fiaGJAM_run3

type <- rep("CA",ncol(y.gjam))

#Running GJAM
rl  <- list(r = 5, N = 20)
ml  <- list(ng = 50000, burnin = 25000, typeNames = type, reductList = rl)
out6<- gjam(form, x.gjam, y.gjam, modelList = ml)
save(out6, file = "~/fiaGJAM_longRS.Rdata")
summary(out6)

plotPars <- list(SMALLPLOTS=F, GRIDPLOTS = T, SAVEPLOT = T)
p <-gjamPlot(out6, plotPars)




plot(rowSums(out6$inputs$y), rowSums(out6$prediction$ypredMu), ylab = "Predicted biomass m2/ha", 
     xlab = "Observed biomass m2/ha",main = "Observed vs Predicted")
abline(0, 1, col = "red", lty = 2)

gjam.lm <- lm(out2$prediction$ypredMu ~ out2$inputs$y)

write.csv(out6$inputs$y, "ydataTru_in_FIA_Trees_longRS.csv")
write.csv(out6$prediction$ypredMu, "ydataPred_in_FIA_Trees_longRS.csv")
