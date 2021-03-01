#########################
#This code was written by Lillian Hislop
#2020.02.26
#WSMDP Carbohydrate Analysis 
#########################



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
library(tidyr)
library(stringr)

getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

source("R/CarbEquationComparison.R")
source("R/PredictionCulling.R")
source("R/MixedEndoEqnValidation.R")
source("R/AgPredOutput.R")


PredictionStarchFile <-paste("Data/RawData/wsmdp2021gl","st","plsall.txt", sep = "")
PredictionSugarFile <-paste("Data/RawData/wsmdp2021gl","su","plsall.txt", sep = "")

CarbPredData <- data.frame(AgPredOutput(PredictionStarchFile,PredictionSugarFile))

CarbPredData <- separate(CarbPredData, Samples, into = c("Year","Row","Rep"), sep = "([W.Y.-])", remove = FALSE)
CarbPredData$Location <- "Wi"
CarbPredData$Location[which(grepl("N",CarbPredData$Year))] = "NY"
CarbPredData$Year<-str_remove(CarbPredData$Year,"N")
CarbPredData <- CarbPredData[-c(1:2),]
tail(CarbPredData)
CarbPredData$Year <- as.numeric(CarbPredData$Year)


summary(CarbPredData)
hist(CarbPredData$Year)
hist(CarbPredData$Glucose)
hist(CarbPredData$WSP)
hist(CarbPredData$Starch)
hist(CarbPredData$Fructose)
hist(CarbPredData$Sucrose)
hist(CarbPredData$Total.Sugar)

