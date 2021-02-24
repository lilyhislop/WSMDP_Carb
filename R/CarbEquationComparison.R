CarbEquationComparison <- function(ActualFile, PredictionStarchFile,PredictionSugarFile, set,scale, modeltype){
  ######Lets make this a function######
  ####I want to read in the data and have it spit out all the visualizations and a doc reporting scores of the equation#####
  #read in Actual, which will be a data frame. What carb values did we actually get
  #read in prediction, which will be a data frame. What carb values did we predict
  #read in modeltype which will be a string. How was this eqn made?
  #establish variable names
  dfLabels <- c("Samples","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
  
  ######process the actual wetlab data######
  #Read in actual data
  validationSamples <- read.csv(file = ActualFile)
  #establish actual data in dataframe
  Actual <- validationSamples[c(1,3:9)]
  # Actual <- data.frame(validationSamples$NIR_ID,validationSamples$Starch, validationSamples$Total_Poly,
  # validationSamples$WSP, validationSamples$Ave.Glucose, validationSamples$Ave.Fructose, 
  # validationSamples$Ave.Sucro, validationSamples$Ave.Total)
  colnames(Actual) <- dfLabels
  
  ######process the predicted wetlab data ######
  #read in the predicted wetlab data 
  starch_pred <- read.delim(file = PredictionStarchFile ,header = TRUE,skip =10 , sep= ",")
  sugar_pred <- read.delim(file = PredictionSugarFile,header = TRUE,skip =10 , sep= ",")
  
  #establish predicted in dataframe 
  wetlabPredicted <- data.frame(starch_pred$Sample.Number,starch_pred$Starch, starch_pred$Total_Poly,
                                starch_pred$WSP, sugar_pred$AveGlucose, 
                                sugar_pred$AveFructose, sugar_pred$AveSucrose, 
                                sugar_pred$AveTotalSugar)
  
  #rename variables
  colnames(wetlabPredicted) <- dfLabels
  
  #the NIR scanned the wetlab samples multiple times. Average all the wetlab preditions. This also orders the samples alphanumerically
  wetlabPredicted_aggregate <- wetlabPredicted %>%
    group_by(Samples) %>%
    summarise(Starch = mean(Starch),Total.Polysaccharides = mean(Total.Polysaccharides), 
              WSP = mean(WSP), Glucose =mean(Glucose), Fructose = mean(Fructose),
              Sucrose = mean(Sucrose), Total.Sugar = mean(Total.Sugar))
  
  
  
  Both <- merge(wetlabPredicted_aggregate,Actual, by = "Samples")
  ValidationDF <- data.frame()
  
  # ####################
  # #so the output from the NIR hass two blank lines, and three lines that have no matching wetlab names, due to over scanning. I need to remove these
  # #remove the blank lines. there were the calibration scans
  # wetlabPredicted_ordered <- wetlabPredicted_aggregate[-(1:2),]
  # 
  # #find the samples with no wetlab corresponding lines
  # diff <- setdiff(wetlabPredicted_ordered$Samples, Actual$Samples)
  # #I want to automatically remove the differences, but i cant figure out how. Remove the manually for the time being
  # which(wetlabPredicted_ordered$Samples == diff[3]) #differences at 8,13,22
  # 
  # #alright, these match now
  # wetlabPredicted_ordered <- wetlabPredicted_ordered[-c(8,13,22),]
  # # setdiff(wetlabPredicted_ordered$Samples, Actual$Samples)
  
  #set up in dataframe to be processed 
  # Prediction <- as.data.frame(wetlabPredicted_ordered)
  ###########################
  Prediction <- as.data.frame(Both)
  
  # #establish csv file for recording equation results
  # carbStatsFileName <- paste("Data/",set,"/wsmdp2021",scale,modeltype,set,"_NIR_Eqn_Prediction_Stats.csv", sep = "")
  # file.remove(carbStatsFileName)
  # print(paste(scale,modeltype,set))
  #plot all of the comparisons of actual vs predicted and recprd
  for(i in 2:8){
    carbCompare  <- lm(Prediction[,i]~Actual[,i])
    carbFileName <- paste("Figures/",set,"/wsmdp2021",scale,modeltype,set,colnames(Actual[i]),"_NIR_Eqn_Prediction_vis.png", sep = "")
    png(carbFileName)
    par(mfrow=c(1,1))
    # print(summary(carbCompare))
    rsqua <- summary(carbCompare)$r.squared
    plot(Prediction[,i]~Actual[,i],
         pch = 16,
         xlab = paste(colnames(Actual[i])," wetlab (%)",sep = ""),
         ylab = paste(colnames(Actual[i])," NIR Prediction (%)",sep = ""),
         main = paste("Actual Vs Predicted ",scale,colnames(Actual[i])," r^2 =",trunc(rsqua*10^3)/10^3," using ", modeltype ,sep = ""))
    abline(coefficients(carbCompare), lwd = 2, lty = 2, col = "red")
    # text(15,max(Prediction[,i])-5,labels = paste("r^2 =",trunc(rsqua*10^3)/10^3))
    
    slope <- trunc(10^3*summary(carbCompare)$coefficients[2])/10^3
    intercept <- trunc(10^3*summary(carbCompare)$coefficients[1])/10^3
    r <- trunc(10^3*summary(carbCompare)$r.squared)/10^3
    trait <- colnames(Actual[i])
    stats <- data.frame(set, scale, modeltype,trait, slope, intercept, r)
    # ValidationDF <- rbind(ValidationDF, stats)
    
    dev.off()
    
    # sink(carbStatsFileName, append = TRUE)
    # print(colnames(Actual[i]))
    # print(summary(carbCompare))
    # print(paste("slope=",slope,"intercept=",intercept))
    # print(paste("r2=",r))
    # sink()
    ValidationDF <- rbind(ValidationDF, stats)
    
  }
  return(ValidationDF)
}
