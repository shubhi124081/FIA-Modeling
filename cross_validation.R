###### FIA validation script

#from NASAgjamRun.R script, call out of bag samples 
#30% of all observations were removed for validation 

# oob <- which(!(seq(1, nrow(x)) %in% sample))
x.oob <- x[oob,]
y.oob <- y[oob,]

effort <- list(columns = 1:ncol(y.oob), values = rep.int(1, nrow(y.oob)))
new <- list(ydata = y.oob, xdata = x.oob, effort = effort, nsim = 500)
predg <- gjamPredict(out6, newdata = new)

Testcolsum <- colSums(y.oob, na.rm = T)

if (length(Testcolsum[Testcolsum == 0]) != 0){
  
  predg2 <- as.data.frame(predg$sdList$yMu)
  
  zeroDF <- as.data.frame(matrix(0,nrow(predg2),length(Testcolsum[Testcolsum == 0])))
  names(zeroDF) <- names(Testcolsum[Testcolsum == 0])
  
  predg3 <- cbind(predg2,zeroDF)
  predg3 <- predg3[,names(y.oob)]
  
  youtPred <- rbind(y.oob,predg3)
  
} else {
  
  youtPred <- rbind(y.oob,predg$sdList$yMu)
  
}

tabCV <- as.data.frame(matrix(NA,nrow(y.oob),2))
names(tabCV) <- c("Pred","Obs")
youtPred <- as.data.frame(predg$sdList$yMu)
tabCV$Pred <- rowSums(youtPred)
tabCV$Obs <- rowSums(y.oob)

m1 <- lm(Obs~Pred, data = tabCV)
summary(m1)


plot(tabCV$Obs, tabCV$Pred, main = "Observed vs Predicted")
abline(0, 1, col = "red", lty =2)

#loop for R2 for species fit
r2 <- matrix(0, nrow = ncol(y.oob), ncol = 1)
for(i in 1:ncol(y.oob)){
  
  lm.tmp <- lm(predg$sdList$yMu[,i] ~ y.oob[,i])
  r2[i, 1] <- summary(lm.tmp)$r.squared
}




write.csv(y.oob, "ydataTru_out_FIA_trees_newRS.csv")
write.csv(youtPred, "ydataPred_out_FIA_trees_newRS.csv")
