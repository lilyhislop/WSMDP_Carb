CarbDataFrameVis <- function(CarbDF, TitleQuality = ""){
  VarietyPos <- which(colnames(CarbDF) == "Variety")
  SamplesPos <- which(colnames(CarbDF) == "Samples")
  StarchPos <- which(colnames(CarbDF) == "Starch")
  SugarPos <- which(colnames(CarbDF) == "Total.Sugar")
  
  #Visualize the Total Traits
  tiff(paste("Figures/WSMDP_Carb_AllTrait_",TitleQuality,"_NIREqnOut_Histogram.tiff",sep=""),  width = 750, height = 500)
  par(mfrow=c(3,3))
  for(i in StarchPos:SugarPos){
    hist(pull(CarbDF[,i]),main = paste("Histogram of",colnames(CarbDF[i])), xlab = colnames(CarbDF[i]))
  }
  dev.off()
  
  
  #endomutants 
  tiff(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Starch_byEndo_Boxplot.tiff",sep=""), width = 500, height = 250)
  dat.m <- reshape2::melt(CarbDF, id.vars = "endo", measure.vars = c('Starch', 'WSP', 'Total.Polysaccharides') )
  
    p <- ggplot(dat.m) +
    geom_boxplot(aes(x = endo, y = value, color = variable))
    print( p + labs(title = "Polysaccharide Content by Endosperm Type") + xlab("Endosperm Mutant") + ylab("Kernel Composition (*%*)"))
  dev.off()
  
  
  tiff(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Sugar_byEndo_Boxplot.tiff",sep=""), width = 500, height = 250)
  dat.m <- reshape2::melt(CarbDF, id.vars = "endo", measure.vars = c('Glucose', 'Fructose', 'Sucrose', 'Total.Sugar') )
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = endo, y = value, color = variable))
  print(p + labs(title = "Sugar Content by Endosperm Type") + xlab("Endosperm Mutant") + ylab("Kernel Composition (*%*)"))
  dev.off()
  
  
 
  pdf(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"SugarandStarch_byEndo_Boxplot.pdf",sep=""))
  # par(mfrow=c(1,2))
  dat.m <- reshape2::melt(CarbDF, id.vars = "endo", measure.vars = c('Glucose', 'Fructose', 'Sucrose', 'Total.Sugar') )
  g <- ggplot(dat.m) +
    geom_boxplot(aes(x = endo, y = value, color = variable)) +
    labs(title = "Sugar Content by Endosperm Type") + xlab("Endosperm Mutant") + ylab(expression(paste("Kernel Composition (",italic("%"),")")))
  
  dat.m <- reshape2::melt(CarbDF, id.vars = "endo", measure.vars = c('Starch', 'WSP', 'Total.Polysaccharides') )
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = endo, y = value, color = variable)) +
    labs(title = "Polysaccharide Content by Endosperm Type") + xlab("Endosperm Mutant") + ylab(expression(paste("Kernel Composition (",italic("%"),")")))
  
  figure <- ggarrange(g,p,
                      labels = c("A", "B"),
                      ncol = 1, nrow = 2)
  print(figure)
  dev.off()
  
  
  ###environment
  tiff(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Starch_byEnvi_Boxplot.tiff",sep=""), width = 750, height = 500)
  dat.m <- reshape2::melt(CarbDF, id.vars = "Envi", measure.vars = c('Starch', 'WSP', 'Total.Polysaccharides') )
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = Envi, y = value, color = variable))
  print( p + labs(title = paste("Polysaccharide Content by Environment",TitleQuality)) + xlab("Enviroment") + ylab("Kernel Composition (*%*)"))
  dev.off()
  
  tiff(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"Sugar_byEnvi_Boxplot.tiff",sep=""), width = 1000, height = 500)
  dat.m <- reshape2::melt(CarbDF, id.vars = "Envi", measure.vars = c('Glucose', 'Fructose', 'Sucrose', 'Total.Sugar') )
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = Envi, y = value, color = variable))
  print(p + labs(title = paste("Sugar Content by Environment",TitleQuality)) + xlab("Enviroment") + ylab("Kernel Composition (*%*)"))
  dev.off()
  
  
  
  pdf(paste("Figures/WSMDP_AllNIRPred_",TitleQuality,"SugarandStarch_byEnvi_Boxplot.pdf",sep=""))
  
  dat.m <- reshape2::melt(CarbDF, id.vars = "Envi", measure.vars = c('Starch', 'WSP', 'Total.Polysaccharides') )
  p <- ggplot(dat.m) +
    geom_boxplot(aes(x = Envi, y = value, color = variable))+ 
    labs(title = "Polysaccharide Content by Environment") + xlab("Enviroment") + ylab(expression(paste("Kernel Composition (",italic("%"),")")))
  
  dat.m <- reshape2::melt(CarbDF, id.vars = "Envi", measure.vars = c('Glucose', 'Fructose', 'Sucrose', 'Total.Sugar') )
  g <- ggplot(dat.m) +
    geom_boxplot(aes(x = Envi, y = value, color = variable))+
    labs(title = "Sugar Content by Environment") + xlab("Enviroment") + ylab(expression(paste("Kernel Composition (",italic("%"),")")))
  
  figure <- ggarrange(g,p,
                      labels = c("A", "B"),
                      ncol = 1, nrow = 2,
                      font.label = list(size = 24))
  print(figure)
  dev.off()
  
}