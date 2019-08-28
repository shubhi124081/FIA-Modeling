library(reshape2)
library(gjam)
library(stringr)
library(dplyr)
library(rgdal)
library(rgeos)


#convert outputs to a grid

#this example uses the 4km grid cell, but the BBS grid cell parameters will be different

library(foreign)
library(rgdal)
library(raster)

pixelSize <- 0.041665991

predx2 <- readOGR(dsn = paste0(getwd(),"/Grid4kmUSAcont/shape"), layer = "Grid4kmUSAcont")

rowRas <- ((max(as.numeric(as.character(predx2@data$LAT)))+(pixelSize/2))-(min(as.numeric(as.character(predx2@data$LAT)))-(pixelSize/2)))/(pixelSize)
colRas <- ((max(as.numeric(as.character(predx2@data$LON)))+(pixelSize/2))-(min(as.numeric(as.character(predx2@data$LON)))-(pixelSize/2)))/(pixelSize)

rasterBlank <- raster(nrow=rowRas, 
                      ncol=colRas, 
                      xmn=min(as.numeric(as.character(predx2@data$LON)))-(pixelSize/2), 
                      xmx=max(as.numeric(as.character(predx2@data$LON)))+(pixelSize/2), 
                      ymn=min(as.numeric(as.character(predx2@data$LAT)))-(pixelSize/2), 
                      ymx=max(as.numeric(as.character(predx2@data$LAT)))+(pixelSize/2), 
                      crs = CRS(proj4string(predx2)))
filesL <- list.files(path ='Z:/R/FIAgjam/Pred')

outRootName <- "final"
sp.info <- read.csv("Z:/R/FIAgjam/REF_SPECIES.csv", header = T)
pb <- txtProgressBar(max = 100)
for (i in 1:1){#
  print("reading in csv...")
  inputTab <- read.csv('Z:/R/FIAgjam/Pred/out4_FIA_hist.csv')
  
  inputTab[,1:(ncol(inputTab)-1)] <- inputTab[,1:(ncol(inputTab)-1)]*(1/4)
  
  print("merging...")
  
  inputTabxy <- merge(predx2,inputTab, by = c("geeID"), all.x = T)
  
  varsT <- names(inputTab)[1:(ncol(inputTab)-1)]
  
  mm <- match(varsT, sp.info$SPECIES_SYMBOL)
  genus <- sp.info$GENUS[mm]
  species <- sp.info$SPECIES[mm]
  
  sp <- data.frame(cbind(varsT, paste0(genus, species)))
  colnames(sp) <- c("varsT", "genusSpecies")
  
  scenario <- 'rcp85'
  timePeriod <- '2040_2069'
  
  print("entering j loop...")
  for (j in 1:length(varsT)){

    outname <- paste0("R:/clark/clark.unix/PBGJAM_Maps/Trees/Climate/PredRas/hist","/R_",outRootName,"_",sp$genusSpecies[j],".tif")
    ras2 <- rasterize(inputTabxy,rasterBlank, field = varsT[j],fun = mean, background = NA, filename = outname, na.rm = T,overwrite=T)
    setTxtProgressBar(pb, value = (j/length(varsT))*100)
  }
}

outpath <- 'R:/clark/clark.unix/PBGJAM_Maps/'
dataTypeMod <- 'Trees/Climate'
scen <- c('rcp45', 'rcp85')
timeP <- c('2040_2069', '2070_2099')
pb <- txtProgressBar(max = 100)

#Create difference maps for future
for (i in 2:length(scen)){
  for (j in 1:length(timeP)){
    
    print(paste("Scenario = ", scen[i]))
    print(paste("Time period = ", timeP[j]))
    
    filesL <- list.files(path = paste0(outpath,dataTypeMod, "/PredRas/",scen[i],"/",timeP[j]),
                         pattern = glob2rx(paste0("R_",outRootName,"_*.tif")))
    
    print("Entering k loop...")
    
    for (k in 1:length(filesL)){
      
      outname2 <- paste0(outpath,dataTypeMod, "/PredRasDif/",scen[i],"/",timeP[j],"/",filesL[k])
      rHist <- raster(paste0(outpath,dataTypeMod,"/PredRas/hist/2018_2018/",filesL[k]))
      rFut <- raster(paste0(outpath,dataTypeMod, "/PredRas/",scen[i],"/",timeP[j],"/",filesL[k]))
      
      rdif <- rFut - rHist
      writeRaster(rdif, outname2, overwrite = T)
      setTxtProgressBar(pb, value = (k/length(filesL))*100)
      
    }    
  }
}

#creates grids for masks

filesL <- list.files(path = paste0(outpath,dataTypeMod,"/PredTab"), pattern = glob2rx(paste0("mask*",dataTypeMod,".csv")))

outRootName <- dataTypeMod

for (i in 1:length(filesL)){#
  inputTab <- read.csv(paste0(outpath,dataTypeMod,"/PredTab/",filesL[i]))
  
  inputTabxy <- merge(predx2,inputTab, by = c("geeID"), all.x = T)
  
  varsT <- names(inputTab)[1:(ncol(inputTab)-1)]
  
  scenario <- unlist(strsplit(filesL[i],"_"))[4]
  timePeriod <- paste0(unlist(strsplit(filesL[i],"_"))[2],"_",unlist(strsplit(filesL[i],"_"))[3])
  
  for (j in 1:length(varsT)){
    outname <- paste0(outpath, dataTypeMod, "/PredRas/",scenario,"/",timePeriod,"/M_",outRootName,"_",varsT[j],".tif")
    ras2 <- rasterize(inputTabxy,rasterBlank, field = varsT[j],fun = mean, background = NA, filename = outname, na.rm = T,overwrite=T)
  }
}

#Cluster analysis
library(vegan)
library(mclust)

#rcp45
inputTab <- read.csv(paste0(outpath,dataTypeMod,"/PredTab/","pred_2018_2018_hist_NEON_Small-Mammals_",dataTypeMod,".csv"))
inputTab$timeID <- "2018_2018_hist"
inputTab$geeIDtimeID <-paste0(inputTab$geeID,"_",inputTab$timeID)

filesL <- list.files(path = paste0(outpath,dataTypeMod,"/PredTab"), pattern = glob2rx(paste0("pred_*rcp45*",dataTypeMod,".csv")))

for (i in 1:length(filesL)){
  inputTab2 <- read.csv(paste0(outpath,dataTypeMod,"/PredTab/",filesL[i]))
  inputTab2$timeID <- paste0(unlist(strsplit(filesL[i],"_"))[2],"_",
                             unlist(strsplit(filesL[i],"_"))[3],"_",
                             unlist(strsplit(filesL[i],"_"))[4])
  inputTab2$geeIDtimeID <-paste0(inputTab2$geeID,"_",inputTab2$timeID)
  inputTab <- rbind(inputTab,inputTab2)
}

##for testing
#sample1 <- sample(1:nrow(inputTab),10000,replace = FALSE)
#inputTabSamp <- inputTab[sample1,]

tabC1 <- inputTab[,1:(ncol(inputTab)-4)]*(1/4)#convert to relative abundance counts/ha
tabC2 <- decostand(tabC1, "norm")
gc()

fit <- mclustBIC(tabC2, G=1:20) #30-40
write.csv(fit[1:50],paste0(outpath,dataTypeMod,"/PredTabCluster/","BIC_1_20_allModels.csv"), row.names = F)

#fit <- mclustBIC(tabC2, G=1:30, modelName = 'VEV') #30-40
#fit <- mclustBIC(tabC2, G=30:45, modelName = 'VEV')

fit <- Mclust(tabC2, G=44, modelName = 'VEV')
inputTab$cluster <- fit$classification
write.csv(inputTab,paste0(outpath,dataTypeMod,"/PredTabCluster/","combo_45_44.csv"), row.names = F)

#summarize mean abundance for each cluster/time period
inputTab$clusterYr <- paste0(inputTab$cluster,".",inputTab$timeID)
inputTabA <- aggregate(.~clusterYr, data = inputTab[,1:(ncol(inputTab)-6)], FUN = mean, na.rm = T)
write.csv(inputTabA,paste0(outpath,dataTypeMod,"/PredTabCluster/","combo_45_44Agg.csv"), row.names = F)

#create rasters of communities

timeP <- unique(inputTab$timeID)

for (k in 1:length(timeP)){
  
  inputTabxy <- merge(predx2,inputTab[inputTab$timeID == timeP[k],c("geeID","cluster")], by = c("geeID"), all.x = T)
  
  scenario <- unlist(strsplit(timeP[k],"_"))[4]
  timePeriod <- paste0(unlist(strsplit(timeP[k],"_"))[2],"_",unlist(strsplit(timeP[k],"_"))[3])
  
  outname <- paste0(outpath, dataTypeMod, "/PredRas/",scenario,"/",timePeriod,"/R_",outRootName,"_cluster_",timePeriod,"_",scenario,".tif")
  ras2 <- rasterize(inputTabxy,rasterBlank, field = "cluster",fun = mean, background = NA, filename = outname, na.rm = T,overwrite=T)
  
}
