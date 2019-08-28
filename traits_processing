#process GJAM results 
betaTable <- out5$parameters$betaTable

xvar <- c("intercept", "prEto", "tmin", "clay", "sand", "soilDepth",
          "u1", "ELEV")

ysp <- colnames(out4$inputs$y)

rfBetaTable <- matrix(NA, nrow = length(ysp), ncol = length(xvar))
colnames(rfBetaTable) <- xvar  
rownames(rfBetaTable) <- ysp
rw = 1
col = 1
for(i in 1:nrow(betaTable)){
  
  #grab covariate
  
  rn <- rownames(betaTable)
  char <- regexpr("_", rn[i])
  xsubstr <- substr(rn[i], char[1] + 1, nchar(rn[i]))
  
  #grab species
  
  ysubstr <- substr(rn[i], 1, char[1] - 1)
  
  #input to table
  if(ysubstr == rownames(rfBetaTable)[rw]){
    print(paste0(ysubstr, ": Correct match for species"))
  
  if(xsubstr == colnames(rfBetaTable)[col]){
    print(paste0(xsubstr, ": Correct match for covariate"))
  
  
  }
  }
  rfBetaTable[rw, col] <- betaTable[i, 1]
  #update column and row
  col = col + 1
  
  if( i %% 8 == 0){
    rw = rw + 1
    col = 1
  }
  print(paste("I is" , i))
  # print(paste("Row is" , rw))
  # print(paste("Col is ", col))
}

#sensitivity table by covariate

sensCov <- function(xvar, ysp, betaTable){
  RFBetaTable <- matrix(NA, nrow = length(ysp), ncol = ncol(betaTable) -1 )
  colnames(RFBetaTable) <- colnames(betaTable)[-5]
  rownames(RFBetaTable) <- ysp
  rw = 1
  i= 1
  
  for(i in 1:nrow(betaTable)){
    
    #grab covariate
    rn <- rownames(betaTable)
    char <- regexpr("_", rn[i])
    xsubstr <- substr(rn[i], char[1] + 1, nchar(rn[i]))
    
    #grab species
    ysubstr <- substr(rn[i], 1, char[1] - 1)
    
    
    #input to table
    if(xsubstr == xvar){
      
      print(paste0(ysubstr, ": Correct match for species"))
    for(j in 1:(ncol(betaTable) - 1)){
        RFBetaTable[rw,j ] <- betaTable[i,j]
      
    }
    
      rw = rw + 1
      }
    # print(rw); print(i)
    }
    
  
  return(as.data.frame(RFBetaTable))
 
  }
  
prSens <- sensCov(xvar = "MOD13Q1_EVI", ysp = ysp, betaTable = betaTable)
tminSens <- sensCov(xvar = "tmin", ysp = ysp, betaTable = betaTable)

#plots
require(ggplot2)

range = 1:30

g <- ggplot(prSens[range,], aes(x = as.factor(rownames(prSens[range,])), y = prSens[range, 1])) + 
  geom_errorbar(aes(ymax = prSens[range,4], ymin = prSens[range, 3])) +
  geom_abline(slope = 0, intercept = 0, colour = "red") +geom_point() + theme_bw()

#colour outliers 

h <- summary(prSens[, 1])

prSens$Outlier <- rep(0, nrow(prSens))
#lower 25% labelled "1"
prSens$Outlier <- ifelse(prSens$Estimate < as.numeric(h[2]), 1, prSens$Outlier) 
#higher 25% labelled "2"
prSens$Outlier <- ifelse(prSens$Estimate > as.numeric(h[5]), 2, prSens$Outlier) 


g <- ggplot(prSens[range,], aes(x = as.factor(rownames(prSens[range,])), y = prSens[range, 1],  
                                colour = as.factor(prSens[range, 5]))) + 
  geom_errorbar(aes(ymax = prSens[range,4], ymin = prSens[range, 3])) +
  geom_abline(slope = 0, intercept = 0, colour = "red") +geom_point() + theme_bw()


#read in trait table 

trait <- read.csv("TraitTable.csv", stringsAsFactors = F)

sp.info$traitNames<- paste0(tolower(sp.info$GENUS), tolower(sp.info$SPECIES))

mm <- match(trait$genusSpecies, sp.info$traitNames,)

spcd <- sp.info$SPECIES_SYMBOL[mm]

trait$SPCD <- spcd

subTrait <- trait[, subsetCols(c("SPCD","maxht", "leafN", "leafP", "SLA", "gmPerSeed", "gmPerCm",
                                 "xylem", "shade", "drought", "flood"), trait)]


prSens$SPCD <- rownames(prSens)


prSens1 <- left_join(prSens, subTrait, by = "SPCD")


prSens_c <- prSens1[which(!(prSens1$Outlier == 0)), ]
range = 1:nrow(prSens_c)

g <- ggplot(prSens_c[range,], aes(x = as.factor((prSens_c$SPCD[range])), y = prSens_c[range, 1],  
                                colour = as.factor(prSens_c$drought[range]))) + 
  geom_errorbar(aes(ymax = prSens_c[range,4], ymin = prSens_c[range, 3])) +
  geom_abline(slope = 0, intercept = 0, colour = "red") +geom_point() + scale_color_hue("Drought") + theme_bw()

#files to write out for NASA


col <- which(names(prSens1) == "SPCD")
traitSp <- prSens1[, col:ncol(prSens1)]
sp.info$SPECIES_SYMBOL <- as.character(sp.info$SPECIES_SYMBOL)

mm <- match(traitSp$SPCD, sp.info$SPECIES_SYMBOL)
traitSp$species <- sp.info$SPECIES[mm]
traitSp$Genus <- sp.info$GENUS[mm]

write.csv(traitSp, "traits_FIA_trees_newRS.csv")

#beta table 

char <- matrix(0, ncol = 1, nrow = length(rownames(betaTable)))
rn <- rownames(betaTable)
for(i in 1:length(rownames(betaTable))){
  
  
  tmp <- regexpr("_", rn[i])
  char[i, 1] <- tmp[1]
}

#grab species 

betaTableOut <- betaTable

betaTableOut$species <- substr(rn, 1, char - 1)
betaTableOut$Xvar <- substr(rn, char + 1, nchar(rn))

rownames(betaTableOut) <- NULL


write.csv(betaTableOut, "betas_FIA_trees_newRS.csv")

