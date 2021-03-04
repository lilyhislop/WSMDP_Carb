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
library(tidyr)
library(stringr)
library(psych)
library(reshape2)#for multiple groups in GGplot


getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

#Read in packages made for these codesets
source("R/CarbEquationComparison.R")
source("R/PredictionCulling.R")
source("R/MixedEndoEqnValidation.R")
source("R/AgPredOutput.R")
source("R/CarbOutlierCleanup.R")


#Read in the output from NIR
NIRStarchOutput <-paste("Data/RawData/wsmdp_allsamples_inclval_wsmdp2021gl","st","pls.txt", sep = "")
NIRSugarOutput <-paste("Data/RawData/wsmdp_allsamples_inclval_wsmdp2021gl","su","pls.txt", sep = "")

#Read in the sample information
SampleInfo <- read.csv("Data/WSMDP_2014-2015_WINY_SampleInfo.csv")
#Fix the sample numbers that are labeled in the NY style e.g " 14A0255" or "15A0005" and change to just the row numbers
SampleInfo$Row <- SampleInfo$Sample.Number
SampleInfo$Row[which(nchar(SampleInfo$Row)>6)] <- substr(SampleInfo$Row[which(nchar(SampleInfo$Row)>6)],4,15)
tail(SampleInfo)
#Generate NIR Code
SampleInfo$NIRBase <- paste(substr(SampleInfo$year,3,4), ifelse(SampleInfo$location == "Wisconsin", "W","NY"), SampleInfo$Row,sep = "")
tail(SampleInfo)

#Assemble that together into a dataframe
CarbNIREqnOut <- data.frame(AgPredOutput(NIRStarchOutput,NIRSugarOutput,condense = TRUE))

#how does that look
head(CarbNIREqnOut)

#Cleanup the NIR output and add rows corresponding to their year and location
CarbNIREqnOut <- separate(CarbNIREqnOut, Samples, into = c("Year","Row","Rep"), sep = "([W.Y.-])", remove = FALSE)
CarbNIREqnOut$Location <- "W"
CarbNIREqnOut$Location[which(grepl("N",CarbNIREqnOut$Year))] = "NY"
CarbNIREqnOut$Year<-str_remove(CarbNIREqnOut$Year,"N")
CarbNIREqnOut <- CarbNIREqnOut[-c(1:2),]
tail(CarbNIREqnOut)
# CarbNIREqnOut$Year <- as.numeric(CarbNIREqnOut$Year)
CarbNIREqnOut$NIRBase <- paste(CarbNIREqnOut$Year,CarbNIREqnOut$Location, CarbNIREqnOut$Row, sep = "")
#Visualize the distribution of the traits from the NIR predicted outputs
summary(CarbNIREqnOut)
colnames(CarbNIREqnOut)

#Combine the Predicted outcome from the NIR with the information about the samples

setdiff(CarbNIREqnOut$NIRBase,SampleInfo$NIRBase)
setdiff(SampleInfo$NIRBase,CarbNIREqnOut$NIRBase)
CarbInfoExpandedDF <- merge(CarbNIREqnOut,SampleInfo, by = "NIRBase")
head(CarbInfoExpandedDF)
summary(CarbInfoExpandedDF)


#Visualize the Total Traits
png("Figures/WSMDP_Carb_AllTrait_NIREqnOut_wOutliers_Histogram.png")
par(mfrow=c(3,3))
for(i in 6:12){
hist(CarbInfoExpandedDF[,i],main = paste("Histogram of",colnames(CarbNIREqnOut[i])), xlab = colnames(CarbNIREqnOut[i]))
}
dev.off()

png("Figures/WSMDP_AllNIRPred_Starch_byEndo_wOUtliers_Boxplot.png", width = 750, height = 500)
dat.m <- melt(CarbInfoExpandedDF, id.vars = "endo", measure.vars = c('Starch', 'WSP', 'Total.Polysaccharides') )
p <- ggplot(dat.m) +
  geom_boxplot(aes(x = endo, y = value, color = variable))
p + labs(title = "Polysaccharide Content by Endosperm Type") + xlab("Endosperm Mutant") + ylab("Percentage of Kernel")
dev.off()


png("Figures/WSMDP_AllNIRPred_Sugar_byEndo_wOutliers_Boxplot.png", width = 1000, height = 500)
dat.m <- melt(CarbInfoExpandedDF, id.vars = "endo", measure.vars = c('Glucose', 'Fructose', 'Sucrose', 'Total.Sugar') )
p <- ggplot(dat.m) +
  geom_boxplot(aes(x = endo, y = value, color = variable))
p + labs(title = "Sugar Content by Endosperm Type") + xlab("Endosperm Mutant") + ylab("Percentage of Kernel")
dev.off()

CleanedInfo <- CarbOutlierCleanup(CarbInfoExpandedDF)



