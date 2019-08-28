

setwd("R:/clark/clark.unix/GEEtables/final/FIA/out")

filesL <- dir(getwd())
tmpkeep <- c( "MOD13Q1_EVI_monthly.csv", "MOD13Q1_NDVI_monthly.csv", 
              "MOD16A2_ET_monthly.csv", "TRMM_prSum_monthly.csv",
             "TRMM_pr_monthly.csv")
mm <- match( tmpkeep, filesL)
filesKeep <- filesL[mm]

# smos <- read.csv(file = paste0(getwd(),"/", "SMOS_smp_monthly.csv"))
# 
# smos$year <- substr(smos$yr_m, 1, 4)
# 
# smos  <- smos[, -which(names(smos) == "yr_m")]
# 
# dir.create("Z:/R/AGGsmos")
# setwd("Z:/R/AGGsmos")
# 
# year <- unique(smos$year)
# outPath <- "Z:/R/AGGsmos"



for(j in 5:length(filesKeep)){
  
  pb <- txtProgressBar(max = 100)
  
  print("Reading csv...")
  
  tab <- read.csv(file = paste0(getwd(), "/", filesKeep[j]))
  tab$year <- substr(tab$yr_m, 1, 4)
  tab <- tab[, -which(names(tab) == "yr_m")]
  
  year <- unique(tab$year)
  char <- regexpr("monthly", filesKeep[j])
  varName <- substr(filesKeep[j], 1, char[1]-2)
  tmpPath <- paste0("Z:/R/AGG", varName)
  dir.create(tmpPath)
  
  print("Entering i loop...")  

for(i in 1:length(year)){

  tmp <- tab[tab$year == year[i],]
  
  varCol <- which(names(tab) == varName)
  agg <- aggregate(tmp[, varCol], by = list(tmp$geeID, tmp$year), FUN =mean, na.rm= T)
  
  colnames(agg) <- c("geeID", "year", varName)
  # if(nrow(agg) == length(agg$geeID)){print("Aggregated correctly!")}
  
  write.csv(agg, file = paste0(tmpPath, "/", "AggData_", year[i], ".csv"),row.names = F)

  setTxtProgressBar(pb, value = (i/length(year))*100)
    
}

agtmp <- dir(tmpPath)

print("Entering m loop...")
collect <- list()
for(m in 1:length(agtmp)){
  collect[[m]] <- read.csv(file = paste0(tmpPath, "/", "AggData_", year[m], ".csv"), stringsAsFactors = F, header = T)
}

tmp.yearly <- do.call(rbind, collect)

print("Writing final csv...")
write.csv(tmp.yearly, file = paste0("z:/R/", varName, "_yearly.csv"), row.names = F)
}
