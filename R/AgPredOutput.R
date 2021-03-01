AgPredOutput <- function(PredictionStarchFile,PredictionSugarFile){
  #this function reads in the WinISI output of predicted starch and sugar data and outputs a dataframe of that data combined and cleaned up
  
  
######process the predicted wetlab data ######
#read in the predicted wetlab data 
dfLabels <- c("Samples","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")

starch_pred <- read.delim(file = PredictionStarchFile ,header = TRUE,skip =10 , sep= ",")
sugar_pred <- read.delim(file = PredictionSugarFile,header = TRUE,skip =10 , sep= ",")
SugStar_Pred <- merge(starch_pred,sugar_pred, by = "Sample.Number")

#establish predicted in dataframe 
wetlabPredicted <- data.frame(SugStar_Pred$Sample.Number,SugStar_Pred$Starch, SugStar_Pred$Total_Poly,
                              SugStar_Pred$WSP, SugStar_Pred$AveGlucose, 
                              SugStar_Pred$AveFructose, SugStar_Pred$AveSucrose, 
                              SugStar_Pred$AveTotalSugar)

#rename variables
colnames(wetlabPredicted) <- dfLabels

#the NIR scanned the wetlab samples multiple times. Average all the wetlab preditions. This also orders the samples alphanumerically
wetlabPredicted_aggregate <- wetlabPredicted %>%
  group_by(Samples) %>%
  summarise(Starch = mean(Starch),Total.Polysaccharides = mean(Total.Polysaccharides), 
            WSP = mean(WSP), Glucose =mean(Glucose), Fructose = mean(Fructose),
            Sucrose = mean(Sucrose), Total.Sugar = mean(Total.Sugar))

return(wetlabPredicted_aggregate)
}