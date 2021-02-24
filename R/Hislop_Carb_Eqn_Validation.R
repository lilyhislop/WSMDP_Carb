#########################
#This code was written by Lillian Hislop
#2020.12.11
#WSMDP Carbohydrate Analysis 
#Working with NIR Prediction Data
#########################



#####Todo::::: Make mixed automated for multiple models



#########################
###Establish Workspace###
#########################
#begin by establishing a new fresh work space
rm(list=ls())
library(dplyr)
library(lme4)
library(ggplot2)
library(ggpmisc)
library(grDevices)

#establish the working directory
getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

source("R/CarbEquationComparison.R")
source("R/PredictionCulling.R")
#########################
###Equation Evaluation###
#########################
#read in the weblab of the validation samples
ActualFile <- "Data/WSMDP_EqnValidation_Wetlab_Data.csv"

#######Visualize the different models#####

models <- c("pls","pcr","modpls")
scaleVerb <- c("Global","su1 calibrated","sh2 calibrated")
scaleAbr <- c("gl","su1","sh2")
setVerb <- "Validation"
setArb <- "val"
ValidationDF <- data.frame()
for(j in 1:3){
for(i in 1:3){
PredictionStarchFile <-paste("Data/RawData/wsmdp2021",scaleAbr[j],"st",models[i],setArb,".txt", sep = "")
PredictionSugarFile <-paste("Data/RawData/wsmdp2021",scaleAbr[j],"su",models[i],setArb,".txt", sep = "")

ValidationDF <- rbind(ValidationDF, CarbEquationComparison(ActualFile,PredictionStarchFile,PredictionSugarFile,setVerb,scaleAbr[j], models[i]))

}}

write.csv(ValidationDF,"Data/WSMDP_Equations_fit_Validation_Evaluation.csv")

#### Down here lets make an equation that combines both calibration types##########
setVerb <- "Validation"
setArb <- "val"
scaleAbr <- c("su1","sh2","ot")
modeltype <- "pcr" 
  
  PredictionSU1StarchFile <- paste("Data/RawData/wsmdp2021",scaleAbr[1],"st",modeltype,setArb,".txt", sep = "")
  PredictionSU1SugarFile <- paste("Data/RawData/wsmdp2021",scaleAbr[1],"su",modeltype,setArb,".txt", sep = "")
  PredictionSH2StarchFile <- paste("Data/RawData/wsmdp2021",scaleAbr[2],"st",modeltype,setArb,".txt", sep = "")
  PredictionSH2SugarFile <- paste("Data/RawData/wsmdp2021",scaleAbr[2],"su",modeltype,setArb,".txt", sep = "")
  PredictionSH2SugarFile <- paste("Data/RawData/wsmdp2021",scaleAbr[2],"su",modeltype,setArb,".txt", sep = "")
  PredictionOtherSugarFile <- paste("Data/RawData/wsmdp2021",scaleAbr[3],"su",modeltype,setArb,".txt", sep = "")
 
  #establish variable names
  dfLabels <- c("Samples","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
  dfLabels2 <- c("Samples","Endo","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
  
  ######process the actual wetlab data######
  #Read in actual data
  validationSamples <- read.csv(file = ActualFile)
  #establish actual data in dataframe
  
  Actual <- data.frame(validationSamples$NIR_ID,validationSamples$Endo, validationSamples$Starch, validationSamples$Total_Poly,
                       validationSamples$WSP, validationSamples$Ave.Glucose, validationSamples$Ave.Fructose, 
                       validationSamples$Ave.Sucro, validationSamples$Ave.Total)
  colnames(Actual) <- dfLabels2
  
  Actualsu1 <- Actual[which(Actual$Endo==1),]
  Actualsh2 <- Actual[which(Actual$Endo==2),]
  Actualother <- Actual[which(Actual$Endo==3),]
  

BothSU1 <- PredictionCulling(Actualsu1,PredictionSU1StarchFile,PredictionSU1SugarFile)
BothSH2 <- PredictionCulling(Actualsh2,PredictionSH2StarchFile,PredictionSH2SugarFile)
BothOther <- PredictionCulling(Actualother,PredictionSU1StarchFile,PredictionOtherSugarFile)

AllTogetherNow <- rbind(BothSU1, BothSH2, BothOther)

  #plot all of the comparisons of actualvs predicted and recprd
  for(i in 2:8){
    carbCompare  <- lm(AllTogetherNow[,i]~AllTogetherNow[,i+8])
    carbFileName <- paste("Figures/WSMDP_",colnames(AllTogetherNow[i]),"_EndoSpecific_NIR_Eqn_Prediction_Validity.png", sep = "")
    png(carbFileName)
    par(mfrow=c(1,1))
    summary(carbCompare)
    rsqua <- summary(carbCompare)$r.squared
    plot(AllTogetherNow[,i]~AllTogetherNow[,i+8],
         pch = 16,
         xlab = paste(colnames(AllTogetherNow[i])," wetlab (%)",sep = ""),
         ylab = paste(colnames(AllTogetherNow[i])," NIR Prediction (%)",sep = ""),
         main = paste("Actual Vs Predicted ",colnames(AllTogetherNow[i])," r^2 =",trunc(rsqua*10^3)/10^3," using EndoSpecific" ,sep = ""))
    abline(coefficients(carbCompare), lwd = 2, lty = 2, col = "red")
    # text(15,max(Prediction[,i])-5,labels = paste("r^2 =",trunc(rsqua*10^3)/10^3))
    
    dev.off()
    
    
    slope <- trunc(10^3*summary(carbCompare)$coefficients[2])/10^3
    intercept <- trunc(10^3*summary(carbCompare)$coefficients[1])/10^3
    r <- trunc(10^3*summary(carbCompare)$r.squared)/10^3
    trait <- colnames(Actual[i])
    set <- setArb
    scale <- "mixed"
    stats <- data.frame(set, scale, modeltype, trait, slope, intercept, r)
    ValidationDF <- rbind(ValidationDF, stats)
  }
  
write.csv(ValidationDF, "Data/WSMDP_Equations_fit_Validation_Evaluation.csv")


