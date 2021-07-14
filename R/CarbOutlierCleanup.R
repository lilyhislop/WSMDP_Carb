CarbOutlierCleanup <- function(CarbDataFrame,DFType = NA, alpha = 0.05){
#lets tease out outliers
  #Establish to remove list
  sink(paste("Data/OutputtedData/EdittedSampleRecords_VerboseChanges_",DFType,".txt",sep = ""), split = TRUE)
  toremove <- c()
  
  #make sure the original positions of the samples are recorded before we go futtzing around with them 
  CarbDataFrame$DFPosition <- rownames(CarbDataFrame)
  VarietyPos <- which(colnames(CarbDataFrame) == "Variety")
  SamplesPos <- which(colnames(CarbDataFrame) == "Samples")
  StarchPos <- which(colnames(CarbDataFrame) == "Starch")
  SugarPos <- which(colnames(CarbDataFrame) == "Total.Sugar")
  
  #establish outlier record dataframe
  ProbSamples <- data.frame()
  
  #lets get just the sample, variety and carb info positions
  Subset <- c(SamplesPos,c(StarchPos:SugarPos), VarietyPos)
  
  #For each outlier to be analyzed, lets record what we are doing to it 
  recordKeeping <- function(Action){
    out <- data.frame(CarbDataFrame$DFPosition[i],CarbDataFrame[PosOutliers[i],Subset],Action)
    names(out) <- names(ProbSamples)
    ProbSamples <- rbind(ProbSamples,out)
    return(ProbSamples)
  }
##### Remove samples with incorrect rep types
BadRepPos <- which(CarbDataFrame$SampleRep > 4 )
print(paste("The Sample Replicate number above 4 mean it was taken after fresh harvest time. Remove all of these sample.", length(BadRepPos), "samples removed."))
if(length(BadRepPos) != 0){
CarbDataFrame <- CarbDataFrame[-BadRepPos,]  }

##########Iterate through Carb lower Outliers First##########
Carbs <- c("Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
for(c in 1:length(Carbs)){

  #Carb being analyzed and the relevant summation carb
  if(c < 4){CarbSumPos  = StarchPos+1}
  if(c > 3 ){CarbSumPos = SugarPos}
CarbDFPos <- c+StarchPos - 1


CarbDataFrame <- CarbDataFrame[order(CarbDataFrame[CarbDFPos]),]
PosOutliers <- which(CarbDataFrame[CarbDFPos] < 0)
PosOutliersVar <- CarbDataFrame[PosOutliers,"Variety"]
head(CarbDataFrame[PosOutliers, Subset])
  

counter <- 1
end <- length(PosOutliers)+1
i <- 1

#lets look at each potential outlier within the context of the the other examples of that variety. 
while(counter < end){
  #reporting out whats going on
  print(paste("The ",colnames(CarbDataFrame[CarbDFPos])," content of sample ", CarbDataFrame$DFPosition[PosOutliers[counter]], " is ",round(CarbDataFrame[PosOutliers[counter],CarbDFPos],2),"%.",sep = ""))
  if(is.na(CarbDataFrame[PosOutliers[counter],CarbDFPos]))
  {
    print(paste("The",colnames(CarbDataFrame[CarbDFPos]),"content of sample", CarbDataFrame$DFPosition[PosOutliers[counter]], "is NA. Skip this sample."))
  }
  if(!is.na(PosOutliersVar[[1]][i])){
  #ok, so what is the outlier in the context of the variety
  outlierincontext <- CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[[1]][i]),]
  #how much of an outlier is it?
  outlierData <- grubbs.test(pull(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[[1]][i]),colnames(CarbDataFrame[CarbDFPos])]))
  
  if(!is.na(outlierData$p.value) && !is.na(pull(CarbDataFrame[PosOutliers[counter],CarbDFPos]))){
    #if the outlier has a pvalue of less than alpha
    if(outlierData$p.value < alpha)
    {
      #if the summation stat is also an outlierr, then add the sample to the to remove list. If the summation isn't an outlier, set the carb and summation carb to NA
      #test summation carb (total poly or total sugar) by seeing if its greater than standard deviations from the mean 
      polySD <- sd(pull(outlierincontext[CarbSumPos]),na.rm = TRUE)
      polymean <- mean(pull(outlierincontext[CarbSumPos]),na.rm = TRUE)
      OutlierPoly <- pull(CarbDataFrame[PosOutliers[counter],CarbSumPos])

      #Check if the carb's summation stat is also an outlier. (Total.Poly or Total.Sugar) It is more than 3 sd from the mean?
      if((polymean+3*polySD < OutlierPoly || polymean-3*polySD > OutlierPoly) && !is.na(OutlierPoly)){
      
        ProbSamples <- recordKeeping("Removed")
        toremove <- c(toremove,CarbDataFrame$DFPosition[i])
        CarbDataFrame <- CarbDataFrame[-PosOutliers[i],]
        print(paste("It is an outlier among all",PosOutliersVar[[1]][i], "lines. It's ",colnames(outlierincontext[CarbSumPos])," content is also an outlier. The Sample will be removed"))
        
        #To account for the sample being removed, i and the end needs to be adjusted
        counter = counter-1
        end = end -1
      }
      
      # The total summationation stat is not an outlier
      if((polymean+3*polySD > OutlierPoly && polymean-3*polySD < OutlierPoly) || is.na(OutlierPoly)){
        ProbSamples <- recordKeeping(paste(colnames(CarbDataFrame[CarbDFPos]),"Nulled"))
        CarbDataFrame[PosOutliers[i],CarbDFPos]<- NA
        CarbDataFrame[PosOutliers[i],CarbSumPos]<- NA
        print(paste("It is an outlier among all",PosOutliersVar[[1]][i], "lines. It's ",colnames(outlierincontext[CarbSumPos])," content is not an outlier. The Sample will have its ",colnames(CarbDataFrame[CarbDFPos])," and ",colnames(outlierincontext[CarbSumPos])," set to NA"))
        }
      
    }
    #if the outlier has a pvalue larger than the alpha, set the the carb to 0
    if(outlierData$p.value > alpha){
      ProbSamples <- recordKeeping(paste(colnames(CarbDataFrame[CarbDFPos]),"Zeroed"))
      CarbDataFrame[PosOutliers[i],CarbDFPos]<-0
      print(paste("It is not an outlier among all",PosOutliersVar[[1]][i], "lines.",colnames(CarbDataFrame[CarbDFPos])," has been set to 0%, since biologically it can't be less than 0."))
      }
  }
  if(is.na(outlierData$p.value)){
    ProbSamples <- recordKeeping(paste(colnames(CarbDataFrame[CarbDFPos]),"Zeroed"))
    CarbDataFrame[PosOutliers[i],CarbDFPos]<-0
    print(paste("There are not enough",PosOutliersVar[[1]][i], "samples to determine if it is an outlier. ",colnames(CarbDataFrame[CarbDFPos])," has been set to 0%, since biologically it can't be less than 0."))
  }
  }
  if(is.na(PosOutliersVar[[1]][i]))
  {   print(paste("The variety of the sample is NA. Skip this sample."))
  }
  i = i+1
  counter = counter + 1
}

}

###########now lets look at the upper outliers##################
Carbs <- c("Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
  for(c in 1:length(Carbs)){
    
    #Carb being analyzed and the relevant summation carb
    if(c < 4){CarbSumPos  = StarchPos+1}
    if(c > 3 ){CarbSumPos = SugarPos}
    CarbDFPos <- c+StarchPos - 1
    
    CarbDataFrame <- CarbDataFrame[order(CarbDataFrame[CarbDFPos], decreasing = TRUE),]
    
    #what entries are above the 2.698S SD?
    CarbSummary <- quantile(CarbDataFrame[CarbDFPos],na.rm = TRUE)
    IQR <- CarbSummary[4]-CarbSummary[2]
    UpperWisk <- CarbSummary[4]+1.5*IQR
    
    #ok, what are the samples that are above that upper wisker
    PosOutliers <- which(CarbDataFrame[CarbDFPos] > UpperWisk)
    PosOutliersVar <- CarbDataFrame[PosOutliers,"Variety"]
    head(CarbDataFrame[PosOutliers, Subset])
    
    counter <- 1
    end <- length(PosOutliers)+1
    i <- 1
    
    #lets look at each potential outlier within the context of the the other examples of that variety. 
    while(counter < end){
      
      #reporting out whats going on
      print(paste("The ",colnames(CarbDataFrame[CarbDFPos])," content of sample ", CarbDataFrame$DFPosition[PosOutliers[counter]], " is ",round(CarbDataFrame[PosOutliers[counter],CarbDFPos],2),"%. The upper wisker is at ", round(UpperWisk,2),"%.",sep = ""))
      if(is.na(CarbDataFrame[PosOutliers[counter],CarbDFPos]))
      {
        print(paste("The",colnames(CarbDataFrame[CarbDFPos]),"content of sample", CarbDataFrame$DFPosition[PosOutliers[counter]], "is NA. Skip this sample."))
      }
      if(!is.na(PosOutliersVar[[1]][i])){
      #ok, so what is the outlier in the context of the variety
      outlierincontext <- CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[[1]][i]),]
      #how much of an outlier is it?
      outlierData <- grubbs.test(pull(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[[1]][i]),colnames(CarbDataFrame[CarbDFPos])]))
      
      if(!is.na(outlierData$p.value ) && !is.na(pull(CarbDataFrame[PosOutliers[counter],CarbDFPos]))){
        #if the outlier has a pvalue of less than alpha
        if(outlierData$p.value < alpha)
        {
          #if the summation stat is also an outlierr, then add the sample to the to remove list. If the summation isn't an outlier, set the carb and summation carb to NA
          #test summation carb (total poly or total sugar) by seeing if its greater than standard deviations from the mean 
          polySD <- sd(pull(outlierincontext[CarbSumPos]),na.rm = TRUE)
          polymean <- mean(pull(outlierincontext[CarbSumPos]),na.rm = TRUE)
          OutlierPoly <- pull(CarbDataFrame[PosOutliers[counter],CarbSumPos])
          
          #Check if the carb's summation stat is also an outlier. (Total.Poly or Total.Sugar) It is more than 3 sd from the mean?
          if((polymean+3*polySD < OutlierPoly || polymean-3*polySD > OutlierPoly) && !is.na(OutlierPoly)){
            
            ProbSamples <- recordKeeping("Removed")
            toremove <- c(toremove,CarbDataFrame$DFPosition[i])
            CarbDataFrame <- CarbDataFrame[-PosOutliers[i],]
            print(paste("It is an outlier among all",PosOutliersVar[[1]][i], "lines. It's ",colnames(outlierincontext[CarbSumPos])," content is also an outlier. The Sample will be removed"))
            
            #To account for the sample being removed, i and the end needs to be adjusted
            counter = counter-1
            end = end -1
          }
          
          # The total summationation stat is not an outlier
          if((polymean+3*polySD > OutlierPoly && polymean-3*polySD < OutlierPoly) || is.na(OutlierPoly)){
            ProbSamples <- recordKeeping(paste(colnames(CarbDataFrame[CarbDFPos]),"Nulled"))
            CarbDataFrame[PosOutliers[i],CarbDFPos]<- NA
            CarbDataFrame[PosOutliers[i],CarbSumPos]<- NA
            print(paste("It is an outlier among all",PosOutliersVar[[1]][i], "lines. It's ",colnames(outlierincontext[CarbSumPos])," content is not an outlier. The Sample will have its ",colnames(CarbDataFrame[CarbDFPos])," and ",colnames(outlierincontext[CarbSumPos])," set to NA."))
          }
          
        }
        #if the outlier has a pvalue larger than the alpha, set the the carb to 0
        if(outlierData$p.value > alpha){
          # ProbSamples <- recordKeeping(paste(colnames(CarbDataFrame[CarbDFPos]),"Zeroed"))
          # CarbDataFrame[PosOutliers[i],CarbDFPos]<-0
          print(paste("It is not an outlier among all",PosOutliersVar[[1]][i], "lines.",colnames(CarbDataFrame[CarbDFPos])," has not been altered."))
        }
      }
      if(is.na(outlierData$p.value)){
        # ProbSamples <- recordKeeping(paste(colnames(CarbDataFrame[CarbDFPos]),"Zeroed"))
        # CarbDataFrame[PosOutliers[i],CarbDFPos]<-0
        print(paste("There are not enough",PosOutliersVar[[1]][i], "samples to determine if it is an outlier. ",colnames(CarbDataFrame[CarbDFPos])," has not been altered."))
      }
      }
      if(is.na(PosOutliersVar[[1]][i]))
      {   print(paste("The variety of the sample is NA. Skip this sample."))
      }
      i = i+1
      counter = counter + 1
    
    }
}
colnames(ProbSamples) <- c("Position", colnames(CarbDataFrame[PosOutliers[i],Subset]),"Edit")
write.csv(file = paste("Data/OutputtedData/EdittedSampleRecords_",DFType,".csv",sep = ""),ProbSamples)
sink()

return(CarbDataFrame)
}
