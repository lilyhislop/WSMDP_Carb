#########################
#This code was written by Lillian Hislop
#2020.12.11
#WSMDP Carbohydrate Analysis 
#Working with NIR Prediction Data
#########################



#####Todo::::: Compare actual vs predicted. where are the predictions of the actual?
#make validation set more robust
#statistically ennumerate how good the equation is
#output the best equation into a CSV


#########################
###Establish Workspace###
#########################
#begin by establishing a new fresh work space
rm(list=ls())
library(dplyr)
# library(lme4)
# library(ggplot2)
# library(ggpmisc)
# library(grDevices)

#establish the working directory
getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

source("R/CarbEquationComparison.R")
source("R/PredictionCulling.R")
source("R/MixedEndoEqnValidation.R")
source("R/AgPredOutput.R")

#########################
###Equation Evaluation###
#########################
#read in the weblab of the validation samples
ActualFile <- "Data/WSMDP_EqnValidation_Wetlab_Data.csv"

#######Visualize the Validation Dataset for the the different models#####
####this code looks at all validation data set for the different models and outputs their R^2 value
#all the differen models
models <- c("pls","pcr","modpls")

#the different calibration sets
scaleVerb <- c("Global","su1 calibrated","sh2 calibrated")
scaleAbr <- c("gl","su1","sh2")

#Val equations
setVerb <- "Validation"
setArb <- "val"
ValidationDF <- data.frame()

#itterate through the different models and calibration types 
for(j in 1:3){
for(i in 1:3){
PredictionStarchFile <-paste("Data/RawData/wsmdp2021",scaleAbr[j],"st",models[i],setArb,".txt", sep = "")
PredictionSugarFile <-paste("Data/RawData/wsmdp2021",scaleAbr[j],"su",models[i],setArb,".txt", sep = "")

#read out the outputs of those comparisons
ValidationDF <- rbind(ValidationDF, CarbEquationComparison(ActualFile,PredictionStarchFile,PredictionSugarFile,setVerb,scaleAbr[j], models[i]))
}}

#what if we only apply the sh2 calibration to the sh2 lines, etc. This code does that
ValidationDF <- MixedEndoEqnValidation(ActualFile,ValidationDF)

#write all of that to a csv
write.csv(ValidationDF, "Data/OutputtedData/WSMDP_Equations_fit_Validation_Evaluation.csv")

#Which model is the best over all?
ValDFAve <- aggregate(r~scale + modeltype, ValidationDF, FUN = mean)
ValDFAve[order(ValDFAve$r,decreasing = TRUE)[1],] #this one is the highest scoring model! 



###########statistical Comparison of predicted vs actual for the samples used to create eqn
ComparisonSu1SugarFile <- read.delim("Data/RawData/wsmdp_allsamples_inclval_wwet_wsmdp2021su1supls.txt", header = TRUE, sep = "", skip = 18)
ComparisonSh2SugarFile <- read.delim("Data/RawData/wsmdp_allsamples_inclval_wwet_wsmdp2021sh2supls.txt", header = TRUE, sep = "", skip = 18)
ComparisonOtSugarFile <- read.delim("Data/RawData/wsmdp_allsamples_inclval_wwet_wsmdp2021otsupls.txt", header = TRUE, sep = "", skip = 18)
ComparisonSu1StarchFile <- read.delim("Data/RawData/wsmdp_allsamples_inclval_wwet_wsmdp2021su1stpls.txt", header = TRUE, sep = "", skip = 18)
ComparisonSh2StarchFile <- read.delim("Data/RawData/wsmdp_allsamples_inclval_wwet_wsmdp2021sh2stpls.txt", header = TRUE, sep = "", skip = 18)

#su1
condenseSug <- function(ComparisonDF){
spacers <- which(ComparisonSu1SugarFile[,2] == "vs.")
ComparisonDFFruct <- ComparisonDF[c(1:spacers[1]),2:4]
ComparisonDFGluc <- ComparisonDF[(spacers[1]+8):(spacers[2]-1),2:4]
ComparisonDFSuc <- ComparisonDF[(spacers[2]+8):(spacers[3]-1),2:4]
ComparisonDFTot <- ComparisonDF[(spacers[3]+8):(dim(ComparisonDF)[1]),2:4]

ComparisonSugar<- merge(merge(merge(ComparisonDFFruct, ComparisonDFGluc, by= "Sample"),ComparisonDFSuc, by = "Sample"), ComparisonDFTot, by= "Sample")
colnames(ComparisonSugar) <- c("Samples", "Fructose.Wetlab","Fructose.Eqn", "Glucose.Wetlab","Glucose.Eqn",  "Sucrose.Wetlab","Sucrose.Eqn",  "TotalSug.Wetlab","TotalSug.Eqn")
ComparisonSugar[2:9] <- sapply(ComparisonSugar[2:9],as.numeric)
return(ComparisonSugar)
}
condenseStarch <- function(ComparisonDF){
  spacers <- which(ComparisonDF[,2] == "vs.")
  ComparisonDFStar <- ComparisonDF[1:(spacers[1]-1),2:4]
  ComparisonDFPoly <- ComparisonDF[(spacers[1]+8):(spacers[2]-1),2:4]
  ComparisonDFWSP <- ComparisonDF[(spacers[2]+8):(dim(ComparisonDF)[1]),2:4]
  
  ComparisonStar<- merge(merge(ComparisonDFStar, ComparisonDFPoly, by= "Sample"),ComparisonDFWSP, by = "Sample")
  colnames(ComparisonStar) <- c("Samples", "Starch.Wetlab","Starch.Eqn", "TotalPoly.Wetlab","TotalPoly.Eqn",  "WSP.Wetlab","WSP.Eqn")
  ComparisonStar[2:7] <- sapply(ComparisonStar[2:7],as.numeric)
  return(ComparisonStar)
}

Su1Sug <- condenseSug(ComparisonSu1SugarFile)
Sh2Sug <- condenseSug(ComparisonSh2SugarFile)
OtSug <- condenseSug(ComparisonOtSugarFile)

Su1St <- condenseStarch(ComparisonSu1StarchFile)
Sh2St <- condenseStarch(ComparisonSh2StarchFile)

#now match up these wetlab and equation outputs with the carb info and mash it all up into one

CarbInfosu1s <- CarbInfoExpandedDF[which(CarbInfoExpandedDF$endo == "su1"),]
CarbInfosh2s <- CarbInfoExpandedDF[which(CarbInfoExpandedDF$endo == "sh2"),]
CarbInfoots <- CarbInfoExpandedDF[which(CarbInfoExpandedDF$endo != "sh2" & CarbInfoExpandedDF$endo != "su1" ),]


Su1s <- merge(merge(CarbInfosu1s,Su1Sug, by = "Samples"),Su1St, by = "Samples")
Sh2s <- merge(merge(CarbInfosh2s,Sh2Sug, by = "Samples"),Sh2St, by = "Samples")
Ots <- merge(merge(CarbInfoots,OtSug, by = "Samples"),Su1St, by = "Samples")
MixedEqnFull <- rbind(Su1s,Sh2s,Ots)



#Now I have a variable that has all the projected values, for the values used to calibrate the equations. WHat are the statistics on that?
#set carb to a number 1:7. carb <- c(Fructose, Glucose, Sucrose, Total Sugar, Starch, Total Polysaccharide, WSP)
EqnStats <- function(DF){
  Out <- data.frame(Carb = c("Fructose", "Glucose", "Sucrose", "Total Sugar", "Starch", "Total Polysaccharide", "WSP"),
                    RMSEP = rep(NA,7), 
                    bias = rep(NA,7),
                    SEE = rep(NA,7),
                    slope = rep(NA,7),
                    intercept = rep(NA,7),
                    R2 = rep(NA,7))
  dfpos <- c(2,4,6,8,10,12,14)
  for(carb in 1:7){
  Out$RMSEP[carb] <- sqrt(sum((DF[,dfpos[carb]]- DF[,dfpos[carb]-1])^2)/dim(DF)[1])
  Out$bias[carb] <- mean(DF[,dfpos[carb]]) - mean(DF[,dfpos[carb]-1])
  Out$SEE[carb] <- sqrt((dim(DF)[1]/(dim(DF)[1]-1))*(Out$RMSEP[carb]^2-Out$bias[carb]^2))
  Out$Carb[carb] <- colnames(DF[dfpos[carb]])
  }
  return(Out)
}


#Visualize these
R2Vis <- function(DF, label, Out){
  dfpos <- c(2,4,6,8,10,12,14)
  Carb = c("Fructose", "Glucose", "Sucrose", "Total Sugar", "Starch", "Total Polysaccharide", "WSP")
for(i in 1:7){
  carbCompare  <- lm(DF[,dfpos[i]]~ DF[,dfpos[i]-1])
  carbFileName <- paste("Figures/wsmdp2021_",label,Carb[i],"_NIR_Eqn_Prediction_vis.png", sep = "")
  png(carbFileName)
  par(mfrow=c(1,1))
  print(summary(carbCompare))
  rsqua <- summary(carbCompare)$r.squared
  plot(DF[,dfpos[i]]~ DF[,dfpos[i]-1],
       pch = 16,
       xlab = paste(Carb[i]," wetlab (%)",sep = ""),
       ylab = paste(Carb[i]," NIR Prediction (%)",sep = ""),
       main = paste("Actual Vs Predicted ",Carb[i]," r^2 =",trunc(rsqua*10^3)/10^3,sep = ""))
  abline(coefficients(carbCompare), lwd = 2, lty = 2, col = "red")
  # text(15,max(Prediction[,i])-5,labels = paste("r^2 =",trunc(rsqua*10^3)/10^3))
  
  Out$slope[i] <- trunc(10^3*summary(carbCompare)$coefficients[2])/10^3
  Out$intercept[i] <- trunc(10^3*summary(carbCompare)$coefficients[1])/10^3
  Out$R2[i] <- trunc(10^3*summary(carbCompare)$r.squared)/10^3

  
  dev.off()
}
  return(Out)
}

MixedEqnStats <- EqnStats(MixedEqnFull[,22:35])
MixedEqnStatsR <- R2Vis(MixedEqnFull[,22:35], "Uncleaned_PredVsWetlab_for_Calibration_Samples", MixedEqnStats)
write.csv(MixedEqnStats, "Data/OutputtedData/WSMDP_MixedEquation_Uncleaned_Statistics.csv")


#the NIR scanned the wetlab samples multiple times. Average all the wetlab preditions. This also orders the samples alphanumerically
MixedEqnFullMerged <- MixedEqnFull %>%
  group_by(Samples) %>%
  summarise(Starch = mean(Starch),Total.Polysaccharides = mean(Total.Polysaccharides), 
            WSP = mean(WSP), Glucose =mean(Glucose), Fructose = mean(Fructose),
            Sucrose = mean(Sucrose), Total.Sugar = mean(Total.Sugar),
            Fructose.Wetlab = mean(Fructose.Wetlab), Fructose.Eqn = mean(Fructose.Eqn),
            Glucose.Wetlab = mean(Glucose.Wetlab), Glucose.Eqn = mean(Glucose.Eqn),
            Sucrose.Wetlab = mean(Sucrose.Wetlab), Sucrose.Eqn = mean(Sucrose.Eqn),
            TotalSug.Wetlab = mean(TotalSug.Wetlab), TotalSug.Eqn = mean(TotalSug.Eqn),
            Starch.Wetlab = mean(Starch.Wetlab), Starch.Eqn = mean(Starch.Eqn),
            TotalPoly.Wetlab = mean(TotalPoly.Wetlab), TotalPoly.Eqn = mean(TotalPoly.Eqn),
            WSP.Wetlab = mean(WSP.Wetlab), WSP.Eqn = mean(WSP.Eqn))

MixedEqnStats <- EqnStats(data.frame(MixedEqnFullMerged[,9:22]))
MixedEqnStatsR <- R2Vis(data.frame(MixedEqnFullMerged[,9:22]), "Merged_PredVsWetlab_for_Calibration_Samples", MixedEqnStats)
write.csv(MixedEqnStats, "Data/OutputtedData/WSMDP_MixedEquation_Merged_Statistics.csv")
