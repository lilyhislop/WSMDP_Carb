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
library(lattice)

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
#6 files types. The starch and sugar predictions and the equations made by calibrating with only high wsp lines, high wsp var sug, high wsp high sug
CarbCombos <- c("hwspsu","hwspst","lwspwfsu","lwspwfst","lwspnfsu","lwspwfst")

#establish DF to hold the file readins 
CarbDF <- list()

#iterate throught the 5 file types and read them in. skip the first 9 lines which are all file descriptors, no the data
for(i in 1:6){CarbDF[[i]] <- read.csv(paste("Data/RawData/wsmdp",CarbCombos[i],"pls_predA-G.txt",sep = ""), skip = 9)}

########now match up these predicted values with the sample info and mash it all up into one########
#combine the starch and sugar files for each equation type. Condense takes all the repeated scans of the same sample and averages them
HWSPsDF <- AgPredOutput(StarchDF = CarbDF[[1]], SugDF = CarbDF[[2]], condense = TRUE)
LWSPWFsDF <- AgPredOutput(StarchDF = CarbDF[[3]], SugDF = CarbDF[[4]], condense = TRUE)
LWSPNFDF <- AgPredOutput(StarchDF = CarbDF[[5]], SugDF = CarbDF[[6]], condense = TRUE)

#Add the endosperm and variety info from the SampleInfo DF to each DF
HWSPsDF$endo <- SampleInfo$endo[match(HWSPsDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
HWSPsDF$Variety <- SampleInfo$Variety[match(HWSPsDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
LWSPWFsDF$endo <- SampleInfo$endo[match(LWSPWFsDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
LWSPWFsDF$Variety <- SampleInfo$Variety[match(LWSPWFsDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
LWSPNFDF$endo <- SampleInfo$endo[match(LWSPNFDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]
LWSPNFDF$Variety <- SampleInfo$Variety[match(LWSPNFDF$NIRBase, SampleInfo$NIRBase, nomatch = NA)]

#Eliminate irrelivant endosperm mutants from each df. We don't want to include sh2 samples that were predicted by the HWSP calibrated equations
#should sh2i samples be estimated by the sh2 calibrated equations or the other calibrated equation?
HWSPs <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se"),]
LWSPWFs <- LWSPWFsDF[which(LWSPNFDF$endo != "se" & LWSPNFDF$endo != "su1"),]
LWSPNFs <- LWSPNFDF[which(LWSPNFDF$endo == "sh2" | LWSPNFDF$endo == "sh2i"),]

#mash it all together 
CarbInfoExpandedWFDF <- rbind(HWSPs,LWSPWFs)
CarbInfoExpandedNFDF <- rbind(HWSPs,LWSPNFs)
CarbInfoExpandedWFDF$Envi <- paste(CarbInfoExpandedWFDF$Year,CarbInfoExpandedWFDF$Location)
CarbInfoExpandedNFDF$Envi <- paste(CarbInfoExpandedNFDF$Year,CarbInfoExpandedNFDF$Location)
CarbDataFrameVis(CarbInfoExpandedWFDF,"WithField_WithOutliers")
CarbDataFrameVis(CarbInfoExpandedNFDF,"NoField_WithOutliers")

CleanedInfoWF <- CarbOutlierCleanup(CarbInfoExpandedWFDF,alpha = 0.05)
CleanedInfoNF <- CarbOutlierCleanup(CarbInfoExpandedNFDF,alpha = 0.05)
CarbDataFrameVis(CleanedInfoWF,"WithField_Cleaned")
CarbDataFrameVis(CleanedInfoNF,"NoField_Cleaned")

#write the names of the varieties used to a csv file so we can find the corresponding GBS data
write.csv(file = "Data/OutputtedData/InbredsWithinWSMDPCarbData.csv",unique(CleanedInfoNF[c("Variety", "endo")]))




##########Linear Model Analysis!##########
head(CleanedInfoWF)

VarDF <- data.frame("Carb" = colnames(CleanedInfoWF)[5:11],"Variety" = rep(NA,7),"Envi" = rep(NA,7), "Rep" = rep(NA,7),"endo" = rep(NA,7),"Variety:Envi" = rep(NA,7),"Residuals"= rep(NA,7))

for(i in 1:7){
formula1 <- paste(colnames(CleanedInfoWF)[i+4],"~ Variety*Envi + Rep + endo")
fit1 <- lm(formula1,data=CleanedInfoWF)
AIC1 <- extractAIC(fit1)
AIC1
summary(fit1)$r.square
print(anova(fit1))
out <- anova(fit1)
SStotal <- sum(out$`Sum Sq`)
for(j in 1:6){
VarDF[i,j+1] <- (out$`Sum Sq`[j]/SStotal)}
}
VarDFMelt <- melt(VarDF)
# # VarDFMelt %>%
# #   gather(Carb,value) %>%
# # ggplot(aes(x = value, y = Carb)) + 
# #   geom_col(position = "identity")#+
# #   #facet_wrap(~variable) #+
# # # coord_flip()
# 
# p <- ggplot(data = VarDF, aes(x = Carb, y = Variety)) +
#   geom_bar(stat = "identity")
# p + coord_flip()
# 
#  ggplot(VarDFMelt, aes(x =Carb, y = value)) +
#   facet_wrap(~variable)+
#   geom_bar(aes(fill = factor(Carb)))
# p + coord_flip()

png(paste("Figures/WSMDP_AllNIRPred_MixedEqn_PercentVarianceExplainedby_Factors_WithField.png",sep=""), width = 1000, height = 500)
barchart(~value|variable, group = factor(Carb), data= VarDFMelt,main = "Percent Phenotypic Variance Explained",layout = c(6,1),
         key = simpleKey(text = colnames(CleanedInfoWF)[5:11],
                         rectangles = TRUE, points = FALSE, space = "right"))
dev.off()


length(which(CleanedInfoWF$endo   == "field"))
