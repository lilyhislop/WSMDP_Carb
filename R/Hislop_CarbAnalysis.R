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
library(ggplot2)
library(ggpmisc)
library(grDevices)
library(tidyr)
library(stringr)
library(psych)

getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

#Read in packages made for these codesets
source("R/CarbEquationComparison.R")
source("R/PredictionCulling.R")
source("R/MixedEndoEqnValidation.R")
source("R/AgPredOutput.R")

#Read in the output from NIR
NIRStarchOutput <-paste("Data/RawData/wsmdp2021gl","st","plsall.txt", sep = "")
NIRSugarOutput <-paste("Data/RawData/wsmdp2021gl","su","plsall.txt", sep = "")

#Read in the sample information
SampleInfo <- read.csv("Data/WSMDP_2014-2015_WINY_SampleInfo.csv")
#Fix the sample numbers that are labeled in the NY style e.g " 14A0255" or "15A0005" and change to just the row numbers
SampleInfo$Row <- SampleInfo$Sample.Number
SampleInfo$Row[which(nchar(SampleInfo$Row)>6)] <- substr(SampleInfo$Row[which(nchar(SampleInfo$Row)>6)],4,15)
tail(SampleInfo)
#Generate NIR Code
SampleInfo$NIRBase <- paste(substr(SampleInfo$year,3,4), ifelse(SampleInfo$location == "Wisconsin", "W","NY"), SampleInfo$Row,sep = "")

#Assemble that together into a dataframe
CarbNIREqnOut <- data.frame(AgPredOutput(NIRStarchOutput,NIRSugarOutput,condense = TRUE))

#how does that look
head(CarbNIREqnOut)

#Cleanup the NIR output and add rows corresponding to their year and location
CarbNIREqnOut <- separate(CarbNIREqnOut, Samples, into = c("Year","Row","Rep"), sep = "([W.Y.-])", remove = FALSE)
CarbNIREqnOut$Location <- "Wi"
CarbNIREqnOut$Location[which(grepl("N",CarbNIREqnOut$Year))] = "NY"
CarbNIREqnOut$Year<-str_remove(CarbNIREqnOut$Year,"N")
CarbNIREqnOut <- CarbNIREqnOut[-c(1:2),]
tail(CarbNIREqnOut)
CarbNIREqnOut$Year <- as.numeric(CarbNIREqnOut$Year)

#Visualize the distribution of the traits from the NIR predicted outputs
summary(CarbNIREqnOut)
colnames(CarbNIREqnOut)
png("Figures/WSMDP_Carb_AllTrait_NIREqnOut_Histogram.png")
par(mfrow=c(3,3))
for(i in 5:11){
hist(CarbNIREqnOut[,i],main = paste("Histogram of",colnames(CarbNIREqnOut[i])), xlab = colnames(CarbNIREqnOut[i]))
}
dev.off()


