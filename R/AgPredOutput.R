AgPredOutput <- function(PredictionStarchFile = NA,PredictionSugarFile = NA,StarchDF=NA, SugDF=NA, condense = TRUE){
  #this function reads in the WinISI output of predicted starch and sugar data and outputs a dataframe of that data combined and cleaned up
  
  
######process the predicted wetlab data ######
#read in the predicted wetlab data 
dfLabels <- c("Samples","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
if(!is.na(PredictionStarchFile)){
starch_pred <- read.delim(file = PredictionStarchFile ,header = TRUE,skip =10 , sep= ",")
sugar_pred <- read.delim(file = PredictionSugarFile,header = TRUE,skip =10 , sep= ",")}

if(is.data.frame(StarchDF)){
  starch_pred <- StarchDF
  sugar_pred <- SugDF
}

SugStar_Pred <- merge(starch_pred,sugar_pred, by = "Sample.Number")

#establish predicted in dataframe 
wetlabPredicted <- data.frame(SugStar_Pred$Sample.Number,SugStar_Pred$Starch, SugStar_Pred$Total_Poly,
                              SugStar_Pred$WSP, SugStar_Pred$AveGlucose, 
                              SugStar_Pred$AveFructose, SugStar_Pred$AveSucrose, 
                              SugStar_Pred$AveTotalSugar)

#rename variables
colnames(wetlabPredicted) <- dfLabels

wetlabPredictedOut <- separate(wetlabPredicted, Samples, into = c("Year","Row","Rep"), sep = "([W.Y.-])", remove = FALSE)
wetlabPredictedOut$Location <- "W"
wetlabPredictedOut$Location[which(grepl("N",wetlabPredictedOut$Year))] = "NY"
wetlabPredictedOut$Year<-str_remove(wetlabPredictedOut$Year,"N")

tail(wetlabPredictedOut)
# CarbNIREqnOut$Year <- as.numeric(CarbNIREqnOut$Year)
wetlabPredictedOut$NIRBase <- paste(wetlabPredictedOut$Year,wetlabPredictedOut$Location, wetlabPredictedOut$Row, sep = "")
#Visualize the distribution of the traits from the NIR predicted outputs
summary(wetlabPredictedOut)
colnames(wetlabPredictedOut)
wetlabPredictedOut$Year <- as.integer(wetlabPredictedOut$Year)

if(condense){
#the NIR scanned the wetlab samples multiple times. Average all the wetlab preditions. This also orders the samples alphanumerically
  # first suddenly stopped working. Make a work around
  wetlabPredictedOutSum <- wetlabPredictedOut %>%
  group_by(Samples) %>%
  # summarise(Year = first(Year), Row = first(Row), Rep = first(Rep), Starch = mean(Starch), Total.Polysaccharides = mean(Total.Polysaccharides), 
  #           WSP = mean(WSP), Glucose =mean(Glucose), Fructose = mean(Fructose),
  #           Sucrose = mean(Sucrose), Total.Sugar = mean(Total.Sugar), Location = first(Location), NIRBase = first(NIRBase))
    summarise( Starch = mean(Starch), Total.Polysaccharides = mean(Total.Polysaccharides), 
              WSP = mean(WSP), Glucose =mean(Glucose), Fructose = mean(Fructose),
              Sucrose = mean(Sucrose), Total.Sugar = mean(Total.Sugar))
  
  wetlabPredictedOutSum <- separate(wetlabPredictedOutSum, Samples, into = c("Year","Row","Rep"), sep = "([W.Y.-])", remove = FALSE)
  wetlabPredictedOutSum$Location <- "W"
  wetlabPredictedOutSum$Location[which(grepl("N",wetlabPredictedOutSum$Year))] = "NY"
  wetlabPredictedOutSum$Year<-str_remove(wetlabPredictedOutSum$Year,"N")
  
  tail(wetlabPredictedOutSum)
  # CarbNIREqnOut$Year <- as.numeric(CarbNIREqnOut$Year)
  wetlabPredictedOutSum$NIRBase <- paste(wetlabPredictedOutSum$Year,wetlabPredictedOutSum$Location, wetlabPredictedOutSum$Row, sep = "")
  #Visualize the distribution of the traits from the NIR predicted outputs
  
    }
if(!condense){
  wetlabPredictedOutSum <- wetlabPredictedOut
}
wetlabPredictedOutSum <- wetlabPredictedOutSum[-c(1:2),]
return(wetlabPredictedOutSum)
}