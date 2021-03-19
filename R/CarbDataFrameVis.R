CarbDataFrameVis <- function(CarbDF, TitleQuality = ""){
  VarietyPos <- which(colnames(CarbDF) == "Variety")
  SamplesPos <- which(colnames(CarbDF) == "Samples")
  StarchPos <- which(colnames(CarbDF) == "Starch")
  SugarPos <- which(colnames(CarbDF) == "Total.Sugar")
  
  #Visualize the Total Traits
  png(paste("Figures/WSMDP_Carb_AllTrait_",TitleQuality,"_NIREqnOut_Histogram.png",sep=""))
  par(mfrow=c(3,3))
  for(i in StarchPos:SugarPos){
    hist(pull(CarbDF[,i]),main = paste("Histogram of",colnames(CarbDF[i])), xlab = colnames(CarbDF[i]))
  }
  dev.off()
  
  
  #endomutants 
  png(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Starch_byEndo_Boxplot.png",sep=""), width = 750, height = 500)
  dat.m <- melt(CarbDF, id.vars = "endo", measure.vars = c('Starch', 'WSP', 'Total.Polysaccharides') )
  
    p <- ggplot(dat.m) +
    geom_boxplot(aes(x = endo, y = value, color = variable))
    print( p + labs(title = paste("Polysaccharide Content by Endosperm Type",TitleQuality)) + xlab("Endosperm Mutant") + ylab("Percentage of Kernel"))
  dev.off()
  
  
  png(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Sugar_byEndo_Boxplot.png",sep=""), width = 1000, height = 500)
  dat.m <- melt(CarbDF, id.vars = "endo", measure.vars = c('Glucose', 'Fructose', 'Sucrose', 'Total.Sugar') )
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = endo, y = value, color = variable))
  print(p + labs(title = paste("Sugar Content by Endosperm Type",TitleQuality)) + xlab("Endosperm Mutant") + ylab("Percentage of Kernel"))
  dev.off()
  
  ###environment
  png(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Starch_byEnvi_Boxplot.png",sep=""), width = 750, height = 500)
  dat.m <- melt(CarbDF, id.vars = "Envi", measure.vars = c('Starch', 'WSP', 'Total.Polysaccharides') )
  
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = Envi, y = value, color = variable))
  print( p + labs(title = paste("Polysaccharide Content by Environment",TitleQuality)) + xlab("Enviroment") + ylab("Percentage of Kernel"))
  dev.off()
  
  png(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Sugar_byEnvi_Boxplot.png",sep=""), width = 1000, height = 500)
  dat.m <- melt(CarbDF, id.vars = "Envi", measure.vars = c('Glucose', 'Fructose', 'Sucrose', 'Total.Sugar') )
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = Envi, y = value, color = variable))
  print(p + labs(title = paste("Sugar Content by Environment",TitleQuality)) + xlab("Enviroment") + ylab("Percentage of Kernel"))
  dev.off()
  
}