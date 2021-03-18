#########################
#This code was written by Lillian Hislop
#2020.02.26
#WSMDP Carbohydrate Analysis 
#########################

####What this code hopes to accomplish:
#exploratory work, look for outlines
#visualized data, scatter plot matrix
#output the line names for all the NIR samples, 
#GET gBS DATA for all the nir sample lines
#Anova for the lines
#Do a PCA for all the samples GBS data 
#Run GWAS


#########################
###Establish Workspace###
#########################
#begin by establishing a new fresh work space
rm(list=ls())
library(dplyr)
library(lme4)
library(ggplot2)#for making pretty graphs
library(ggpmisc)
library(grDevices)
library(tidyr) #for separate
library(stringr)#str_remove
library(psych)
library(reshape2)#for multiple groups in GGplot
library(outliers)

getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

#Read in packages made for these codesets
source("R/CarbEquationComparison.R")
source("R/PredictionCulling.R")
source("R/MixedEndoEqnValidation.R")
source("R/AgPredOutput.R")
source("R/CarbOutlierCleanup.R")
source("R/CarbDataFrameVis.R")


#Read in the sample information.
#this file containes the year, location,variety and endosperm of every sample
SampleInfo <- read.csv("Data/WSMDP_2014-2015_WINY_SampleInfo.csv")
#Fix the sample numbers that are labeled in the NY style e.g " 14A0255" or "15A0005" and change to just the row numbers
SampleInfo$Row <- SampleInfo$Sample.Number
SampleInfo$Row[which(nchar(SampleInfo$Row)>6)] <- substr(SampleInfo$Row[which(nchar(SampleInfo$Row)>6)],4,15)
tail(SampleInfo)
#Generate NIR Code. NIR codes are all YearLocationRow - Rep
SampleInfo$NIRBase <- paste(substr(SampleInfo$year,3,4), ifelse(SampleInfo$location == "Wisconsin", "W","NY"), SampleInfo$Row,sep = "")
tail(SampleInfo)


#########Read in the output from NIR##########
####Copied pasted and modified from Hislop_carb_eqn_validation.R
#5 files types. The starch and sugar predictions and the equations made by calibrating with only su1, sh2, and all other sample types
CarbCombos <- c("su1st","sh2st","su1su","sh2su","otsu")

#establish DF to hold the file readins 
CarbDF <- list()

#iterate throught the 5 file types and read them in. skip the first 9 lines which are all file descriptors, no the data
for(i in 1:5){CarbDF[[i]] <- read.csv(paste("Data/RawData/wsmdp_allsamples_inclval_wsmdp2021",CarbCombos[i],"pls.txt",sep = ""), skip = 9)}

########now match up these predicted values with the sample info and mash it all up into one########
#combine the starch and sugar files for each equation type. Condense takes all the repeated scans of the same sample and averages them
Su1sDF <- AgPredOutput(StarchDF = CarbDF[[1]], SugDF = CarbDF[[3]], condense = TRUE)
Sh2sDF <- AgPredOutput(StarchDF = CarbDF[[2]], SugDF = CarbDF[[4]], condense = TRUE)
OtDF <- AgPredOutput(StarchDF = CarbDF[[1]], SugDF = CarbDF[[5]], condense = TRUE)

#Add the endosperm and variety info from the SampleInfo DF to each DF
Su1sDF$endo <- SampleInfo$endo[match(Su1sDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
Su1sDF$Variety <- SampleInfo$Variety[match(Su1sDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
Sh2sDF$endo <- SampleInfo$endo[match(Sh2sDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
Sh2sDF$Variety <- SampleInfo$Variety[match(Sh2sDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
OtDF$endo <- SampleInfo$endo[match(OtDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
OtDF$Variety <- SampleInfo$Variety[match(OtDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]

#Eliminate irrelivant endosperm mutants from each df. We don't want to include sh2 samples that were predicted by the su1 calibrated equations
#should sh2i samples be estimated by the sh2 calibrated equations or the other calibrated equation?
Su1s <- Su1sDF[which(Su1sDF$endo == "su1"),]
Sh2s <- Sh2sDF[which(Sh2sDF$endo == "sh2"),]
Ots <- OtDF[which(OtDF$endo != "sh2" & OtDF$endo != "su1"),]

#mash it all together 
CarbInfoExpandedDF <- rbind(Su1s,Sh2s,Ots)
CarbDataFrameVis(CarbInfoExpandedDF,"WithOutliers")

CleanedInfo <- CarbOutlierCleanup(CarbInfoExpandedDF,alpha = 0.05)
CarbDataFrameVis(CleanedInfo,"Cleaned")




