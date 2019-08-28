#loop to merge all states -----
#FL throwing up errors (10)
#OR throwing up errors (38)
require(dplyr)
#combine all recent trees into one file
setwd('/Users/shubhi/Documents/R/Clark/')
recent.path<- '/Users/shubhi/Documents/R/Clark/FIAMerged3/'
all.recent<- dir(recent.path)



st.names <- substr(all.recent, 0,2)
seq <- seq(0, length(st.names), by = 2)
st.names <- st.names[ seq]
skip.st<- c('AS','GU','FM', 'HI','MP','PR','PW','VI')

all.plots<- all.recent[which(all.recent%in%paste0(st.names, '_allPlots.csv'))]
all.plots <- all.plots[which(!(all.plots%in%paste0(skip.st, '_allPlots.csv')))]

all.tree <- all.recent[which(all.recent%in%paste0(st.names, '_allTrees.csv'))]
all.tree <- all.tree[which(!(all.tree%in%paste0(skip.st, '_allTrees.csv')))]

out.path <- paste0(getwd(), '/FIA_final2/')
#for allPlot and allTree files

for(f in 1:length(all.plots)){
  
  
  #plot files 
  
  tmp.file<- all.plots[f]
  tmp<- read.csv(paste0(recent.path, tmp.file), stringsAsFactors = FALSE, as.is = TRUE)
  
  # last.char <- nchar(names(tmp[ncol(tmp)]))
  # t.no <- as.numeric(substr(names(tmp[ncol(tmp)]),last.char, last.char))
  
  s.no.plt <- unique(substr(names(tmp), nchar(names(tmp)) -1, nchar(names(tmp))))
  s.no1.plt <- paste0("idreal_plt_", s.no.plt)
  start.plt.df <- as.data.frame(matrix(tmp$CN_t1, nrow(tmp), length(s.no.plt)))
  
  
  names(start.plt.df) <- s.no1.plt
  
  total.plot <- cbind(tmp, start.plt.df)
  
  total.time.plot <- list()
  plot.subset <- list()
  
  # plot.keep <- c('CN', 'INVYR', 'STATECD', 'PLOT', 'PLOT_STATUS_CD', 'DESIGNCD',
  #                'LAT', 'LON', 'ELEV', 'STDAGE', 'ASPECT', 'SLOPE', 'PCTFOREST', 'idreal_plt_')
  # w.cols <- names(total.plot[which(grepl(names(total.plot), plot.keep))])
  
  
  subsetCols <- c('CN', 'INVYR', 'STATECD', 'PLOT', 'DESIGNCD', 'LAT', 'LON', 'ELEV', 
                'STDAGE', 'SLOPE', 'ASPECT', 'idreal_plt', 'PCTFOREST', 'PHYSCLCD', 'MANUAL', 'KINDCD')
  
   for(i in 1:length(s.no.plt)){
    keepCols <- names(total.plot)[which((grepl(s.no.plt[i], names(total.plot))))]
    plot.test <- total.plot[, keepCols]
    names(plot.test) <- (substr(names(plot.test), 0, nchar(names(plot.test)) - 3))
    total.time.plot[[i]] <- plot.test[!(is.na(plot.test[,'CN'])),]
    tmp <- total.time.plot[[i]]
    plot.subset[[i]] <- tmp[, subsetCols]
  }
  # Column `FORTYPCD` can't be converted from character to integer
  # Column `ECO_UNIT_PNW` can't be converted from integer to character

  # tmp <- total.time.plot[[1]]
  # total.time.plot[[1]] <- tmp[, -47]
  
   # tmp <- total.time.plot[[4]]
   # total.time.plot[[4]] <- tmp [, -75]
  
 
  plot.final <- do.call(rbind,plot.subset)
  print(paste("Plot iteration= ", f, "; State = ", substr(tmp.file,1,2)))
  
  #tree files 
  
  tree.tmp.file <- all.tree[f]
  tree.tmp <- read.csv(paste0(recent.path, tree.tmp.file), stringsAsFactors = FALSE)
  
  
  s.no <- unique(substr(names(tree.tmp), nchar(names(tree.tmp)) -1, nchar(names(tree.tmp))))
  s.no1 <- paste0("idreal_tree_", s.no)
  start.tree.df <- as.data.frame(matrix(tree.tmp$CN_t1, nrow(tree.tmp), length(s.no)))
  
  names(start.tree.df) <- s.no1
  
  total.tree <- cbind(tree.tmp, start.tree.df)
  total.time.tree <- list()
  tree.subset <- list()
  
  subsetTreeCols <- c('PLT_CN', 'TREE','SUBP','DIST', 'CONDID','STATUSCD','CR','SPCD','AZIMUTH', 'DIA', 'HT', 'TPA_UNADJ', 
                'idreal_tree', 'CCLCD', 'RECONCILECD')
  
  for(i in 1:length(s.no)){
    keepCols <- names(total.tree)[which((grepl(s.no[i], names(total.tree))))]
    tree.test <- total.tree[, keepCols]
    names(tree.test) <- (substr(names(tree.test), 0, nchar(names(tree.test)) - 3))
    total.time.tree[[i]] <- tree.test[!(is.na(tree.test[,'CN'])),]
    tmp <- total.time.tree[[i]]
    tree.subset[[i]] <- tmp[, subsetTreeCols]
  }
  
  
  tree.final <-bind_rows(tree.subset)
  print(paste("Tree iteration= ", f, "; State = ", substr(tmp.file,1,2)))
  
  #merge based on plt_cn use left_join
  start.plt.distinct<- plot.final %>% distinct(CN, .keep_all=T)
  tree.merged<- left_join(tree.final, start.plt.distinct, by=c('PLT_CN' = 'CN'))
  
  
  #subset columns wanted 
 
  
  
  # tree.merged <- tree.merged[, keepCols]
  
  write.csv(tree.merged, file=paste0(out.path, substr(tmp.file,1,2), '_TreePlotMerged.csv'), row.names=F)
  
}

input.path <- '/Users/shubhi/Documents/R/Clark/FIA_final2/'

merged.files <- dir(input.path)
all.files <- list()

for(f in 1:length(merged.files)){
  tmp<- read.csv(paste0(input.path, merged.files[f]),header=T, stringsAsFactors = F)
  tmp.state<- substr(merged.files[f],1,2)
  tmp$STATE= tmp.state
  # tmp<- tmp[which(tmp$CONDPROP_UNADJ==1),]
  all.files[[f]]<- tmp
  print(f)
}

allFIA<- do.call(rbind,all.files)

sp.info<- read.csv('REF_SPECIES.csv',header=T)
mm<- match(allFIA$SPCD,sp.info$SPCD)
genus<- sp.info$GENUS[mm]
species<- sp.info$SPECIES[mm]


symbol <- sp.info$SPECIES_SYMBOL[mm]

allFIA$genus<- genus
allFIA$species<- species
allFIA$SPsymbol <- symbol

n <- names(allFIA)
n[which(n == "idreal_tree")] <- "IDtree"
n[which(n == "idreal_plt")] <- "IDplot"

allFIA1 <- allFIA[, -(which(n == "IDplot"))]
names(allFIA) <- n

allFIA$ASPECT <- round(allFIA$ASPECT, digits = 1)

write.csv(allFIA, "FIA_061219.csv")

unique.cood <- unique(paste0(allFIA$LAT, allFIA$LON))

write.csv(unique.cood, "Unique_cood.csv")
