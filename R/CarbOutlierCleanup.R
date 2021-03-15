CarbOutlierCleanup <- function(CarbDataFrame,CarbColumns = c(6:12) ,alpha = 0.05){

#lets tease out outliers
  #Establish to remove list
  toremove <- c()
  
  VarietyPos <- which(colnames(CarbDataFrame) == "Variety")
  SamplesPos <- which(colnames(CarbDataFrame) == "Samples")
  StarchPos <- which(colnames(CarbDataFrame) == "Starch")
  SugarPos <- which(colnames(CarbDataFrame) == "Total.Sugar")
  
  #establish outlier record dataframe
  ProbSamples <- data.frame()
  # ProbSamples <- data.frame()
  Subset <- c(SamplesPos,c(StarchPos:SugarPos), VarietyPos)
  
  recordKeeping <- function(Action){
    out <- data.frame(CarbDataFrame$DFPosition[i],CarbDataFrame[PosOutliers[i],Subset],Action)
    names(out) <- names(ProbSamples)
    ProbSamples <- rbind(ProbSamples,out)
    return(ProbSamples)
  }
##########Sucrose Outliers First##########
Carbs <- c("Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
for(c in 1:length(Carbs)){

CarbDataFrame$DFPosition <- rownames(CarbDataFrame)
CarbDataFrame <- CarbDataFrame[order(CarbDataFrame$Sucrose),]
PosOutliers <- which(CarbDataFrame$Sucrose < 0)
PosOutliersVar <- CarbDataFrame[PosOutliers,"Variety"]
head(CarbDataFrame[PosOutliers, Subset])
  
#lets look at each potential outlier within the context of the the other examples of that variety. Start with Check1 and Check 2 since those will have the most other samples
CarbDataFrame[which(CarbDataFrame$Variety == "Check2"),Subset] #The -5 sucrose looks like an outlier
grubbs.test(pull(CarbDataFrame[which(CarbDataFrame$Variety == "Check2"),'Sucrose']),opposite = FALSE) 

counter <- 1
end <- length(PosOutliers)+1
i <- 1
#Do the same for the other samples 
while(counter < end){
  #reporting out whats going on
  print(paste("The Sucrose content of sample ", CarbDataFrame$DFPosition[PosOutliers[counter]], " is ",round(CarbDataFrame[PosOutliers[counter],10],2),"%. I is ",counter,sep = ""))
  
  #ok, so what is the outlier in the context of the variety
  outlierincontext <- CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[[1]][i]),Subset]
  #how much of an outlier is it?
  outlierData <- grubbs.test(pull(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[[1]][i]),'Sucrose']))
  
  if(!is.na(outlierData$p.value)){
    #if the outlier has a pvalue of less than alpha
    if(outlierData$p.value < alpha)
    {
      # print(paste(" an outlier among all",PosOutliersVar[i], "lines"))
      #then test the total polysaccharide data
      # outlierPolyData <- grubbs.test(outlierincontext$Total.Polysaccharides)
      
      #if the total polysaccharide data is an outlier, then add the sample to the to remove list. If the poly isn't an outlier, set the sucrose and total sugar to NA
      #test total Polysaccharide by seeing if its greater than standard deviations from the mean 
      polySD <- sd(outlierincontext$Total.Polysaccharides)
      polymean <- mean(outlierincontext$Total.Polysaccharides)
      OutlierPoly <- CarbDataFrame$Total.Polysaccharides[PosOutliers[counter]]
      if(polymean+3*polySD < OutlierPoly || polymean-3*polySD > OutlierPoly){
      
      # if(outlierPolyData$p.value < alpha){
        ProbSamples <- recordKeeping("Removed")
        toremove <- c(toremove,CarbDataFrame$DFPosition[i])
        CarbDataFrame <- CarbDataFrame[-PosOutliers[i],]
        print(paste("It is an outlier among all",PosOutliersVar[[1]][i], "lines. It's Polysaccharide content is also an outlier. The Sample will be removed"))
        #To account for the sample being removed, i and the end needs to be adjusted
        counter = counter-1
        end = end -1
        }
      if(polymean+3*polySD > OutlierPoly && polymean-3*polySD < OutlierPoly){
        ProbSamples <- recordKeeping("Sucrose Nulled")
        CarbDataFrame[PosOutliers[i],c(10:11)]<- NA
        print(paste("It is an outlier among all",PosOutliersVar[[1]][i], "lines. It's Polysaccharide content is not an outlier. The Sample will have its sucrose and total sugar set to NA"))
        }
      
    }
    #if the outlier has a pvalue larger than the alpha, set the sucrose to 0
    if(outlierData$p.value > alpha){
      ProbSamples <- recordKeeping("Sucrose Zeroed")
      CarbDataFrame[PosOutliers[i],10]<-0
      print(paste("It is not an outlier among all",PosOutliersVar[[1]][i], "lines. Sucrose has been set to 0%, since biologically it can't be less than 0."))
      }
  }
  if(is.na(outlierData$p.value)){
    ProbSamples <- recordKeeping("Sucrose Zeroed")
    CarbDataFrame[PosOutliers[i],10]<-0
    print(paste("There are not enough",PosOutliersVar[[1]][i], "samples to determine if it is an outlier. Sucrose has been set to 0%, since biologically it can't be less than 0."))
  }
  i = i+1
  counter = counter + 1
}

##########WSP Outliers First##########
CarbDataFrame <- CarbDataFrame[order(CarbDataFrame$WSP),]
PosOutliers <- which(CarbDataFrame$WSP < 0)
PosOutliersVar <- CarbDataFrame[which(CarbDataFrame$WSP < 0),"Variety"]
head(CarbDataFrame[PosOutliers, Subset])
#Do the same for the other samples 
counter <- 1
end <- length(PosOutliers)+1
i <- 1
#Do the same for the other samples 
while(counter < end){
  #reporting out whats going on
  print(paste("The WSP content of sample ", CarbDataFrame$DFPosition[PosOutliers[counter]], " is ",round(CarbDataFrame[PosOutliers[counter],8],2),"%. I is ",counter,sep = ""))
  
  #ok, so what is the outlier in the context of the variety
  outlierincontext <- CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),Subset]
  #how much of an outlier is it?
  outlierData <- grubbs.test(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),'WSP'])
  
  if(!is.na(outlierData$p.value)){
    #if the outlier has a pvalue of less than alpha
    if(outlierData$p.value < alpha)
    {
      # print(paste(" an outlier among all",PosOutliersVar[i], "lines"))
      #then test the total Glucose data
      # outlierPolyData <- grubbs.test(outlierincontext$Total.Glucoses)
      
      #if the total Glucose data is an outlier, then add the sample to the to remove list. If the poly isn't an outlier, set the WSP and total sugar to NA
      #test total Glucose by seeing if its greater than standard deviations from the mean 
      polySD <- sd(outlierincontext$Glucose)
      polymean <- mean(outlierincontext$Glucose)
      OutlierPoly <- CarbDataFrame$Glucose[PosOutliers[counter]]
      if(polymean+3*polySD < OutlierPoly || polymean-3*polySD > OutlierPoly){
        
        # if(outlierPolyData$p.value < alpha){
        ProbSamples <- recordKeeping("Removed")
        toremove <- c(toremove,CarbDataFrame$DFPosition[i])
        CarbDataFrame <- CarbDataFrame[-PosOutliers[i],]
        print(paste("It is an outlier among all",PosOutliersVar[i], "lines. It's Glucose content is also an outlier. The Sample will be removed"))
        #To account for the sample being removed, i and the end needs to be adjusted
        counter = counter-1
        end = end -1
      }
      if(polymean+3*polySD > OutlierPoly && polymean-3*polySD < OutlierPoly){
        ProbSamples <- recordKeeping("WSP Nulled")
        CarbDataFrame[PosOutliers[i],c(7:8)]<- NA
        print(paste("It is an outlier among all",PosOutliersVar[i], "lines. It's Glucose content is not an outlier. The Sample will have its WSP and total Polysaccharides set to NA"))
      }
      
    }
    #if the outlier has a pvalue larger than the alpha, set the WSP to 0
    if(outlierData$p.value > alpha){
      ProbSamples <- recordKeeping("WSP Zeroed")
      CarbDataFrame[PosOutliers[i],8]<-0
      print(paste("It is not an outlier among all",PosOutliersVar[i], "lines. WSP has been set to 0%, since biologically it can't be less than 0."))
    }
  }
  if(is.na(outlierData$p.value)){
    ProbSamples <- recordKeeping("WSP Zeroed")
    CarbDataFrame[PosOutliers[i],8]<-0
    print(paste("There are not enough",PosOutliersVar[i], "samples to determine if it is an outlier. WSP has been set to 0%, since biologically it can't be less than 0."))
  }
  i = i+1
  counter = counter + 1
}



colnames(ProbSamples) <- c("Position", "NIRBase", "Samples","Starch","Total.Polysaccharides","WSP","Glucose","Fructose","WSP","Total.Sugar","Variety","Edit")
write.csv(file = "Data/OutputtedData/EdittedSampleRecords.csv",ProbSamples)
return(CarbDataFrame)
}