CarbOutlierCleanup <- function(CarbDataFrame,CarbColumns = c(6:12),VarietyPos = 15 ,alpha = 0.05){
#lets tease out outliers
  #Establish to remove list
  toremove <- c()
  
  #establish outlier record dataframe
  ProbSamples <- data.frame()
  # ProbSamples <- data.frame()
  
  recordKeeping <- function(Action){
    out <- data.frame(CarbDataFrame$DFPosition[i],CarbDataFrame[PosOutliers[i],c(1:2,CarbColumns, VarietyPos)],Action)
    names(out) <- names(ProbSamples)
    ProbSamples <- rbind(ProbSamples,out)
    return(ProbSamples)
  }
  ##########Sucrose Outliers First##########
  CarbDataFrame$DFPosition <- rownames(CarbDataFrame)
  CarbDataFrame <- CarbDataFrame[order(CarbDataFrame$Sucrose),]
  PosOutliers <- which(CarbDataFrame$Sucrose < 0)
  PosOutliersVar <- CarbDataFrame[which(CarbDataFrame$Sucrose < 0),"Variety"]

#lets look at each potential outlier within the context of the the other examples of that variety. Start with Check1 and Check 2 since those will have the most other samples
CarbDataFrame[which(CarbDataFrame$Variety == "Check1"),c(1:2,CarbColumns, VarietyPos)] #the -11 value is definitely an outlier
grubbs.test(CarbDataFrame[which(CarbDataFrame$Variety == "Check1"),'Sucrose'])#Grubs test confirms
CarbDataFrame[which(CarbDataFrame$Variety == "Check2"),c(1:2,CarbColumns, VarietyPos)] #The -5 sucrose looks like an outlier
grubbs.test(CarbDataFrame[which(CarbDataFrame$Variety == "Check2"),'Sucrose']) #grub test confirms, less confidently. But other values of entry 113 are also wonked so remove



#Do the same for the other samples 
for(i in 1:length(PosOutliers)){
  #reportingr(i in 1:l out whats going on
  print(paste("The Sucrose content of sample ",CarbDataFrame$DFPosition[i], " is ",round(CarbDataFrame[PosOutliers[i],11],2),"%. ",sep = ""))
  
  #ok, so what is the outlier in the context of the variety
  outlierincontext <- CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),c(1:2,CarbColumns, VarietyPos)]
  #how much of an outlier is it?
  outlierData <- grubbs.test(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),'Sucrose'])
  
  if(!is.na(outlierData$p.value)){
    #if the outlier has a pvalue of less than alpha
    if(outlierData$p.value < alpha)
    {
      # print(paste(" an outlier among all",PosOutliersVar[i], "lines"))
      #then test the total polysaccharide data
      outlierPolyData <- grubbs.test(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),'Total.Polysaccharides'])
      #if the total polysaccharide data is an outlier, then add the sample to the to remove list. If the poly isn't an outlier, set the sucrose and total sugar to NA
      if(outlierPolyData$p.value < alpha){
        ProbSamples <- recordKeeping("Removed")
        toremove <- c(toremove,CarbDataFrame$DFPosition[i])
        CarbDataFrame <- CarbDataFrame[-PosOutliers[i],]
        print(paste("It is an outlier among all",PosOutliersVar[i], "lines. It's Polysaccharide content is also an outlier. The Sample will be removed"))
        }
      if(outlierPolyData$p.value > alpha){
        ProbSamples <- recordKeeping("Sucrose Nulled")
        CarbDataFrame[PosOutliers[i],c(11:12)]<- NA
        print(paste("It is an outlier among all",PosOutliersVar[i], "lines. It's Polysaccharide content is not an outlier. The Sample will have its sucrose and total sugar set to NA"))
        }
      
    }
    #if the outlier has a pvalue larger than the alpha, set the sucrose to 0
    if(outlierData$p.value > alpha){
      ProbSamples <- recordKeeping("Sucrose Zeroed")
      CarbDataFrame[PosOutliers[i],11]<-0
      print(paste("It is not an outlier among all",PosOutliersVar[i], "lines. Sucrose has been set to 0%, since biologically it can't be less than 0."))

      }
  }
  if(is.na(outlierData$p.value)){
    ProbSamples <- recordKeeping("Sucrose Nulled")
    CarbDataFrame[PosOutliers[i],11]<-0
    print(paste("There are not enough",PosOutliersVar[i], "samples to determine if it is an outlier. Sucrose has been set to 0%, since biologically it can't be less than 0."))
    }
}

##########WSP Outliers First##########
CarbDataFrame <- CarbDataFrame[order(CarbDataFrame$WSP),]
PosOutliers <- which(CarbDataFrame$WSP < 0)
PosOutliersVar <- CarbDataFrame[which(CarbDataFrame$WSP < 0),"Variety"]

#Do the same for the other samples 
for(i in 1:length(PosOutliers)){
  print(paste("The WSP content of sample ",CarbDataFrame$DFPosition[i], " is ",round(CarbDataFrame[PosOutliers[i],8],2),"%. ",sep = ""))
  
  #ok, so what is the outlier in the context of the variety
  CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),c(1:2,CarbColumns, VarietyPos)]
  #how much of an outlier is it?
  outlierData <- grubbs.test(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),'WSP'])
  
  if(!is.na(outlierData$p.value)){
    #if the outlier has a pvalue of less than alpha
    if(outlierData$p.value < alpha)
    {
      # print(paste(" an outlier among all",PosOutliersVar[i], "lines"))
      #then test the total polysaccharide data
      outlierPolyData <- grubbs.test(CarbDataFrame[which(CarbDataFrame$Variety == PosOutliersVar[i]),'Glucose'])
      #if the total polysaccharide data is an outlier, then add the sample to the to remove list. If the poly isn't an outlier, set the WSP and total sugar to NA
      if(outlierPolyData$p.value < alpha){
        ProbSamples <- recordKeeping("Removed")
        toremove <- c(toremove,PosOutliers[i])
        CarbDataFrame <- CarbDataFrame[-PosOutliers[i],]
        print(paste("It is an outlier among all",PosOutliersVar[i], "lines. It's Glucose content is also an outlier. The Sample will be removed"))
      }
      if(outlierPolyData$p.value > alpha){
        ProbSamples <- recordKeeping("WSP Nulled")
        CarbDataFrame[PosOutliers[i],c(8:9)]<- NA
        print(paste("It is an outlier among all",PosOutliersVar[i], "lines. It's Glucose content is not an outlier. The Sample will have its WSP and Total Polysaccharide set to NA"))
      }
      
    }
    #if the outlier has a pvalue larger than the alpha, set the WSP to 0
    if(outlierData$p.value > alpha){
      ProbSamples <- recordKeeping("WSP Zeroed")
      CarbDataFrame[PosOutliers[i],9]<-0
      print(paste("It is not an outlier among all",PosOutliersVar[i], "lines. WSP has been set to 0%, since biologically it can't be less than 0."))
      
    }
  }
  if(is.na(outlierData$p.value)){
    ProbSamples <- recordKeeping("WSP Nulled")
    CarbDataFrame[PosOutliers[i],9]<-0
    print(paste("There are not enough",PosOutliersVar[i], "samples to determine if it is an outlier. WSP has been set to 0%, since biologically it can't be less than 0."))
  }
}


colnames(ProbSamples) <- c("Position", "NIRBase", "Samples","Starch","Total.Polysaccharides","WSP","Glucose","Fructose","WSP","Total.Sugar","Variety","Edit")
write.csv(file = "Data/OutputtedData/EdittedSampleRecords.csv",ProbSamples)



return(CarbDataFrame)
}