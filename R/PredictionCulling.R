PredictionCulling <- function(Actual, PredictionStarchFile,PredictionSugarFile){
  ######process the predicted wetlab data ######
  #read in the predicted wetlab data 
  starch_pred <- read.delim(file = PredictionStarchFile ,header = TRUE,skip =10 , sep= ",")
  sugar_pred <- read.delim(file = PredictionSugarFile,header = TRUE,skip =10 , sep= ",")
  
  #establish predicted in dataframe 
  wetlabPredicted <- data.frame(starch_pred$Sample.Number,starch_pred$Starch, starch_pred$Total_Poly,
                                starch_pred$WSP, sugar_pred$AveGlucose, 
                                sugar_pred$AveFructose, sugar_pred$AveSucrose, 
                                sugar_pred$AveTotalSugar)
  #rename variables'
  dfLabels2 <- c("Samples","Endo","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
  dfLabels <- c("Samples","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
  colnames(wetlabPredicted) <-dfLabels
  
  #the NIR scanned the wetlab samples multiple times. Average all the wetlab preditions. This also orders the samples alphanumerically
  wetlabPredicted_aggregate <- wetlabPredicted %>%
    group_by(Samples) %>%
    summarise(Starch = mean(Starch),Total.Polysaccharides = mean(Total.Polysaccharides), 
              WSP = mean(WSP), Glucose =mean(Glucose), Fructose = mean(Fructose),
              Sucrose = mean(Sucrose), Total.Sugar = mean(Total.Sugar))
  colnames(Actual) <- dfLabels2
  colnames(wetlabPredicted_aggregate) <- dfLabels
  
  Both <- merge(wetlabPredicted_aggregate,Actual, by = "Samples")
  
  return(Both)
}