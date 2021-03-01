#########################
#This code was written by Lillian Hislop
#2020.12.11
#WSMDP Carbohydrate Analysis 
#Working with NIR Prediction Data
#########################



#####Todo::::: Compare actual vs predicted. where are the predictions of the actual?



#########################
###Establish Workspace###
#########################
#begin by establishing a new fresh work space
rm(list=ls())
# library(dplyr)
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

ValidationDF <- MixedEndoEqnValidation(ActualFile,ValidationDF)

write.csv(ValidationDF, "Data/WSMDP_Equations_fit_Validation_Evaluation.csv")

#Look at quality of equation
#start with global pls equation

setArb = "all"
i = 1
j = 1
PredictionStarchFile <-paste("Data/RawData/wsmdp2021",scaleAbr[j],"st",models[i],setArb,".txt", sep = "")
PredictionSugarFile <-paste("Data/RawData/wsmdp2021",scaleAbr[j],"su",models[i],setArb,".txt", sep = "")

EqnPre <- data.frame(AgPredOutput(PredictionStarchFile,PredictionSugarFile))

# ActualFile <- "Data/WSMDP_Wetlab_StarchSugarData_FormatedForWinISI_WithR.csv"
# 
# dfLabels <- c("Samples","Starch","Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose","Total.Sugar")
# WetlabSamples <- read.csv(file = ActualFile)
# 
# WetLabStarchFirst <- WetlabSamples[c(1,6:8,2:5)]
# colnames(WetLabStarchFirst) <- dfLabels
# 
# Both <- merge(EqnPre,WetLabStarchFirst, by = "Samples")
