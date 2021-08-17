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
library(PopGenome)
library(data.table) #for fread
# BiocManager::install("VariantAnnotation")
library("VariantAnnotation")
library("snpStats")
library("compiler") #needed to make GAPIT work
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
library("devtools")#for intsalling from github
install_github("jendelman/GWASpoly")
library(GWASpoly)#for running Gwas
# library("bigmemory") #to make a matrix big
# library(rrBLUP)
library(emmeans)
library(bigmemory)
library(MuMIn)
library(SNPRelate)


getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

#Read in packages made for these codesets
source("R/CarbEquationComparison.R")
source("R/PredictionCulling.R")
source("R/MixedEndoEqnValidation.R")
source("R/AgPredOutput.R")
source("R/CarbOutlierCleanup.R")
source("R/CarbDataFrameVis.R")
source("R/GWASPolyVis.R")
source("R/GWASPolyRunner.R")
source("R/WritePhenoGenotoFile.R")
source("R/PCAFigureCreation.R")
source("R/hmpToNumeric.R")

#########################
###Read in sample data###
#########################

#Read in the sample information.
#this file contains information about the superblock, row, column data for each sample
BookInfo <- read.csv("Data/20142015_WSMDP_PlantingInfo_Condensation.csv")
#There are some duplicates in planting info condensation. Remove duplicates
BookInfo <- BookInfo[!duplicated(BookInfo$SampleLabel),]
length(unique(BookInfo$Inbred))

#this file containes the year, location,variety and endosperm of every sample
SampleInfo <- read.csv("Data/WSMDP_2014-2015_WINY_SampleInfo.csv")
#Fix the sample numbers that are labeled in the NY style e.g " 14A0255" or "15A0005" and change to just the row numbers
SampleInfo$Row <- SampleInfo$Sample.Number
SampleInfo$Row[which(nchar(SampleInfo$Row)>6)] <- substr(SampleInfo$Row[which(nchar(SampleInfo$Row)>6)],4,15)
tail(SampleInfo)
#Generate NIR Code. NIR codes are all YearLocationRow - Rep
SampleInfo$NIRBase <- paste(substr(SampleInfo$year,3,4), ifelse(SampleInfo$location == "Wisconsin", "W","NY"), SampleInfo$Row,sep = "")
tail(SampleInfo)
str(unique(SampleInfo$Variety))
summary.factor(SampleInfo$endo)
length(unique(SampleInfo$Variety))
setdiff(unique(BookInfo$Inbred),unique(SampleInfo$Variety))

#########Read in the output from NIR##########
####Copied pasted and modified from Hislop_carb_eqn_validation.R
#6 files types. The starch and sugar predictions and the equations made by calibrating with only high wsp lines, high wsp var sug, high wsp high sug
CarbCombos <- c("hwspsu","hwspst","lwspwfsu","lwspwfst","lwspnfsu","lwspwfst")

#establish DF to hold the file readins 
CarbDF <- list()

#iterate throught the 5 file types and read them in. skip the first 9 lines which are all file descriptors, no the data
for(i in 1:6){CarbDF[[i]] <- read.csv(paste("Data/RawData/wsmdp",CarbCombos[i],"pls_predA-G.txt",sep = ""), skip = 9)}


#########################
###Manipulate and clean sample data###
#########################

########now match up these predicted values with the sample info and mash it all up into one########
#combine the starch and sugar files for each equation type. Condense takes all the repeated scans of the same sample and averages them
HWSPsDF <- AgPredOutput(StarchDF = CarbDF[[1]], SugDF = CarbDF[[2]], condense = TRUE)
LWSPWFsDF <- AgPredOutput(StarchDF = CarbDF[[3]], SugDF = CarbDF[[4]], condense = TRUE)
LWSPNFsDF <- AgPredOutput(StarchDF = CarbDF[[5]], SugDF = CarbDF[[6]], condense = TRUE)

InfoCombination <- function(NIRDF, SampleInfoDF, BookInfoDF){
  OutDF <- NIRDF
  OutDF$endo <- SampleInfoDF$endo[match(NIRDF$NIRBase, SampleInfoDF$NIRBase, nomatch = NA)]
  OutDF$Variety <- SampleInfoDF$Variety[match(NIRDF$NIRBase, SampleInfoDF$NIRBase, nomatch = NA)]
  OutDF$BookInbred <- BookInfoDF$Inbred[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$PlotNum <- NIRDF$Row 
  OutDF$SampleRep <- NIRDF$Rep 
  
  OutDF$superblock <- BookInfoDF$superblock[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$block <- BookInfoDF$block[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$EnviRep <- BookInfoDF$rep[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$Row <- BookInfoDF$Row[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$Col <- BookInfoDF$Col[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$Check <- BookInfoDF$Check[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$IsFill <- BookInfoDF$IsFill[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  OutDF$IsExperimental <- BookInfoDF$IsExperimental[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)]
  #add in a year/location environment factor
  OutDF$Envi <- paste0(BookInfoDF$Year[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)], BookInfoDF$Location[match(NIRDF$NIRBase, BookInfoDF$NIRBase, nomatch = NA)])

  return(OutDF)
}

HWSPsDF <- InfoCombination(HWSPsDF,SampleInfo,BookInfo)
LWSPWFsDF <- InfoCombination(LWSPWFsDF,SampleInfo,BookInfo)
LWSPNFsDF <- InfoCombination(LWSPNFsDF,SampleInfo,BookInfo)

#Eliminate irrelivant endosperm mutants from each df. We don't want to include sh2 samples that were predicted by the HWSP calibrated equations
#should sh2i samples be estimated by the sh2 calibrated equations or the other calibrated equation?
HWSPs <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se" | (is.na(HWSPsDF$endo)& (HWSPsDF$Year== "15"|HWSPsDF$Year == "14"))),]
# HWSPs <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se" | is.na(HWSPsDF$endo)),]
LWSPWFs <- LWSPWFsDF[which(LWSPNFsDF$endo != "se" & LWSPNFsDF$endo != "su1"),]
LWSPNFs <- LWSPNFsDF[which(LWSPNFsDF$endo == "sh2" | LWSPNFsDF$endo == "sh2i"),]

# #lets look at the validation data that I predicted from kahtleen and Jareds data and output it for further analysis
# jared <- HWSPsDF[which(HWSPsDF$Year == "sc"),]
# write.csv(jared, file = "Data/OutputtedData/JaredsPredictedLineswHWSPeqn.csv")
# kathleen <- HWSPsDF[which(HWSPsDF$Year != "13"&HWSPsDF$Year != "14"&HWSPsDF$Year != "15"&HWSPsDF$Year != "sc"),]
# write.csv(kathleen,file = "Data/OutputtedData/KathleenPredictedLineswHWSPeqn.csv")
# 
# ########mash it all together ######
# #All the predictions with the field lines
# CarbInfoExpandedWFDF <- rbind(HWSPs,LWSPWFs)
# #all the predictions without the field lines
# CarbInfoExpandedNFDF <- rbind(HWSPs,LWSPNFs)
# 
# 
# 
# #visualize these data sets pre cleaning
# CarbDataFrameVis(CarbInfoExpandedWFDF,"WithField_WithOutliers")
# CarbDataFrameVis(CarbInfoExpandedNFDF,"NoField_WithOutliers")
# 
# #clean up the predictive data frames. Reassign or delete outliers 
# CleanedInfoWF_wexcess <- CarbOutlierCleanup(CarbInfoExpandedWFDF,"WF",alpha = 0.05)
# CleanedInfoNF_wexcess <- CarbOutlierCleanup(CarbInfoExpandedNFDF,"NF",alpha = 0.05)
# CleanedInfoWF <- subset(CleanedInfoWF_wexcess, !is.na(IsExperimental))
# CleanedInfoNF <- subset(CleanedInfoNF_wexcess, !is.na(IsExperimental))
# write.csv(file = "Data/OutputtedData/CleanedInfoWFOutput.csv",CleanedInfoWF)
# write.csv(file = "Data/OutputtedData/CleanedInfoNFOutput.csv",CleanedInfoNF)
# 
# #revisualize the dataframes 
# CarbDataFrameVis(CleanedInfoWF,"WithField_Cleaned")
# CarbDataFrameVis(CleanedInfoNF,"NoField_Cleaned")
# 
# #write the names of the varieties used to a csv file so we can find the corresponding GBS data
# write.csv(file = "Data/OutputtedData/InbredsWithinWSMDPCarbDataNF.csv",unique(CleanedInfoNF[c("Variety", "endo")]))
# write.csv(file = "Data/OutputtedData/InbredsWithinWSMDPCarbDataWF.csv",unique(CleanedInfoWF[c("Variety", "endo")]))
# 
# #########################
# ###Validate that the NIR Equation is good###
# #########################
# 
# #######Equation Validation!##############
# #Now I have a variable that has all the projected values, for the values used to calibrate the equations. WHat are the statistics on that?
# #set carb to a number 1:7. carb <- c(Fructose, Glucose, Sucrose, Total Sugar, Starch, Total Polysaccharide, WSP)
# #this functions needs to have the data frame set up so the carbs being compared are directly next to eachother 
# EqnStats <- function(DF){
#   
#   #establish dataframe used to record stats
#   Out <- data.frame(Carb = c("Fructose", "Glucose", "Sucrose", "Total Sugar", "Starch", "Total Polysaccharide", "WSP"),
#                     RMSEP = rep(NA,7), 
#                     bias = rep(NA,7),
#                     SEE = rep(NA,7),
#                     slope = rep(NA,7),
#                     intercept = rep(NA,7),
#                     R2 = rep(NA,7))
#   dfpos <- c(2,4,6,8,10,12,14)
#   for(carb in 1:7){
#     #Calculated the RMSEP
#     Out$RMSEP[carb] <- sqrt(sum((DF[,dfpos[carb]]- DF[,dfpos[carb]-1])^2, na.rm = TRUE)/dim(DF)[1])
#     #Calculate the Bias
#     Out$bias[carb] <- mean(DF[,dfpos[carb]], na.rm = TRUE) - mean(DF[,dfpos[carb]-1],na.rm = TRUE)
#     #Calculated the SEE
#     Out$SEE[carb] <- sqrt((dim(DF)[1]/(dim(DF)[1]-1))*(Out$RMSEP[carb]^2-Out$bias[carb]^2))
#     Out$Carb[carb] <- colnames(DF[dfpos[carb]])
#   }
#   return(Out)
# }
# 
# 
# #Visualize these
# R2Vis <- function(DF, label, Out){
#   dfpos <- c(2,4,6,8,10,12,14)
#   Carb = c("Fructose", "Glucose", "Sucrose", "Total Sugar", "Starch", "Total Polysaccharide", "WSP")
#   for(i in 1:7){
#     carbCompare  <- lm(DF[,dfpos[i]]~ DF[,dfpos[i]-1])
#     carbFileName <- paste("Figures/wsmdp2021_",label,Carb[i],"_NIR_Eqn_Prediction_vis.png", sep = "")
#     png(carbFileName)
#     par(mfrow=c(1,1))
#     print(summary(carbCompare))
#     rsqua <- summary(carbCompare)$r.squared
#     plot(DF[,dfpos[i]]~ DF[,dfpos[i]-1],
#          pch = 16,
#          xlab = paste(Carb[i]," wetlab (%)",sep = ""),
#          ylab = paste(Carb[i]," NIR Prediction (%)",sep = ""),
#          main = paste("Actual Vs Predicted ",Carb[i]," r^2 =",trunc(rsqua*10^3)/10^3,sep = ""))
#     abline(coefficients(carbCompare), lwd = 2, lty = 2, col = "red")
#     # text(15,max(Prediction[,i])-5,labels = paste("r^2 =",trunc(rsqua*10^3)/10^3))
#     
#     Out$slope[i] <- trunc(10^3*summary(carbCompare)$coefficients[2])/10^3
#     Out$intercept[i] <- trunc(10^3*summary(carbCompare)$coefficients[1])/10^3
#     Out$R2[i] <- trunc(10^3*summary(carbCompare)$r.squared)/10^3
#     
#     
#     dev.off()
#   }
#   return(Out)
# }
# 
# 
# ######Visualize the validation. Look at all the samples that were wetlabbed and used to calibrate the equation
# wetlab <- read.csv("Data/WSMDP_Wetlab_StarchSugarData_FormatedForWinISI_WithR.csv")
# wetlabDF <- wetlab[,c(1,6:8,3,2,4,5)]
# 
# WLLabels <- c("Samples","Starch_WL","Total.Polysaccharides_WL", "WSP_WL","Glucose_WL","Fructose_WL","Sucrose_WL","Total.Sugar_WL")
# colnames(wetlabDF) <- WLLabels 
# 
# WFWL <- merge(CleanedInfoWF, wetlabDF, by = "Samples")
# NFWL <- merge(CleanedInfoNF, wetlabDF, by = "Samples")
# alternatingorder <- c(1,5,28,6,29,7,30,8,31,9,32,10,33,11,34)
# WFWLdf <- WFWL[,alternatingorder]
# NFWLdf <- NFWL[,alternatingorder]
# WFWLdfEqnStatsR <- R2Vis(WFWLdf[,2:15], "CleanedWSPeqnWF_PredVsWetlab_for_Calibration_Samples", EqnStats(WFWLdf[,2:15]))
# NFWLdfEqnStatsR <- R2Vis(NFWLdf[,2:15], "CleanedWSPeqnNF_PredVsWetlab_for_Calibration_Samples", EqnStats(NFWLdf[,2:15]))
# 
# write.csv(WFWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsAllWetlabWF.csv")
# write.csv(NFWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsAllWetlabNF.csv")
# 
# ######Visualize the validation. Only the samples that were NOT used to create the equation
# valwetlab <- read.csv("Data/WSMDP_EqnValidation_Wetlab_Data.csv")
# valwetlabDF <- valwetlab[-2]
# colnames(valwetlabDF)[1] <- "Samples"
# 
# WFValWL <- merge(CleanedInfoWF, valwetlabDF, by = "Samples")
# NFValWL <- merge(CleanedInfoNF, valwetlabDF, by = "Samples")
# WFValWLdf <- WFValWL[,alternatingorder]
# NFValWLdf <- NFValWL[,alternatingorder]
# WFValWLdfEqnStatsR <- R2Vis(WFValWLdf[,2:15], "CleanedWSPeqnWF_PredVsWetlab_ValidationSubset", EqnStats(WFValWLdf[,2:15]))
# NFValWLdfEqnStatsR <- R2Vis(NFValWLdf[,2:15], "CleanedWSPeqnNF_PredVsWetlab_ValidationSubset", EqnStats(NFValWLdf[,2:15]))
# 
# write.csv(WFValWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationWF.csv")
# write.csv(NFValWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationNF.csv")
# 
# 
# ####With Jared WetlabDATA
# # HWSPsJ <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se" | is.na(HWSPsDF$endo)),]
# # alternatingorder2 <- c(1,5,26,6,27,7,28,8,29,9,30,10,31,11,32)
# # CarbInfoExpandedWFJDF <- rbind(HWSPsJ,LWSPWFs)
# # CarbInfoExpandedNFJDF <- rbind(HWSPsJ,LWSPNFs)
# # WFValWLJ <- merge(CarbInfoExpandedWFJDF, valwetlabDF, by = "Samples")
# # NFValWLJ <- merge(CarbInfoExpandedNFJDF, valwetlabDF, by = "Samples")
# # WFValWLJdf <- WFValWLJ[,alternatingorder2]
# # NFValWLJdf <- NFValWLJ[,alternatingorder2]
# # WFValWLdfJEqnStatsR <- R2Vis(WFValWLJdf[,2:15], "UnCleanedWSPeqnWF_PredVsWetlab_ValidationSubset_WJared", EqnStats(WFValWLdf[,2:15]))
# # NFValWLdfJEqnStatsR <- R2Vis(NFValWLJdf[,2:15], "UnCleanedWSPeqnNF_PredVsWetlab_ValidationSubset_WJared", EqnStats(NFValWLdf[,2:15]))
# # 
# # write.csv(WFValWLdfJEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationWF_WJared.csv")
# # write.csv(NFValWLdfJEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationNF_WJared.csv")
# 
# #########################
# ###Linear Model Analysis###
# #########################
# 
# linearmodel <- function(SampleDFtoModel,TitleAddendum){
#   SampleDFtoModel$superblock <- as.factor(SampleDFtoModel$superblock)
#   SampleDFtoModel$Col <- as.factor(SampleDFtoModel$Col)
#   SampleDFtoModel$Row <- as.factor(SampleDFtoModel$Row)
#   SampleDFtoModel$Year <- as.factor(SampleDFtoModel$Year)
#   SampleDFtoModel$Envi <- as.factor(SampleDFtoModel$Envi)
#   SampleDFtoModel$Check <- as.factor(SampleDFtoModel$Check)
#   SampleDFtoModel$block <- as.factor(SampleDFtoModel$block)
#   SampleDFtoModel$Rep <- as.factor(SampleDFtoModel$Rep)
#   SampleDFtoModel$endo <- as.factor(SampleDFtoModel$endo)
#   SampleDFtoModel$PlotNum <- as.factor(SampleDFtoModel$PlotNum)
#   
#   statsfile <- paste0("Data/OutputtedData/WSMDP_CarbPheno_stats_",TitleAddendum,".txt")
#   cat(paste0("Phenotypic Statistics for the Carbohydrates of ",TitleAddendum," Through Mixed Linear Modeling."), file=statsfile, sep="\n", append=FALSE)
#   
#   blupHolder <- data.frame("Variety" = sort(unique(SampleDFtoModel$BookInbred)))
#   c1 <- which(colnames(SampleDFtoModel)=="Starch")
#   c2 <- which(colnames(SampleDFtoModel)=="Total.Sugar")
#   carbs <- colnames(SampleDFtoModel)[c1:c2]
#   for(j in 1:length(carbs)){
#     
#   #output the statistics about this trait
#     summary <- summary(SampleDFtoModel[,c1+j-1], )
#     out <- capture.output(summary)
#     cat(out, file=statsfile, sep="\n", append=TRUE)
#     
#     #####Modeling
#     modelpastecheck<-  paste0(carbs[j], " ~ (1|BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")
#     
#     # Use check if there enough checks in the model set
#     if(length(unique(SampleDFtoModel$Check))>1){
#     modelpastecheck<-  paste0(carbs[j], " ~ Check + (1|BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")
#     }
#     cat(modelpastecheck, file=statsfile, sep="\n", append=TRUE)
#     
#     
#     model <- lmer(modelpastecheck,
#                   data=SampleDFtoModel, REML = TRUE)
#     # Decreasing stopping tolerances
#     strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
#     if (all(model@optinfo$optimizer=="nloptwrap")) {
#       model <- update(model, control=strict_tol)
#     }
#   RandomEffects <- ranef(model)
#   tempBlup <- data.frame("Variety" = rownames(RandomEffects$BookInbred), "BLUP" = RandomEffects$BookInbred)
#   blupHolder <- merge(blupHolder, tempBlup, by = "Variety", all = TRUE)
#   colnames(blupHolder)[j+1] <- paste0(carbs[j],".BLUP")
#   
#   # CheckAssumptions
#   png(paste0("Figures/WSMDP_LinearModel_assumptions_",TitleAddendum,"_", carbs[j], ".png"), width = 1000, height = 500)
#   par(mfrow=c(1,3))
#   
#   # Model Fit with REML
#   plot(fitted(model), residuals(model), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
#   abline(h=0,col="red", lwd=1, lty=1)
#   # histogram of residuals
#   hist(residuals(model),main=paste0("Histogram of ", carbs[j]," residuals from ", TitleAddendum),freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
#   x=seq(-5e-15,9e-15,5e-15)
#   curve(dnorm(x,mean(residuals(model)),sd(residuals(model))),add=T,lwd=2, col="red", lty=1)
#   # qq plot
#   qqnorm(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles") 
#   
#   dev.off()
#    # Summary of random effects
#   summary <- summary(model, correlation=FALSE)
#   out <- capture.output(summary)
#   cat(out, file=statsfile, sep="\n", append=TRUE)
#   
#   # Write out residuals from ANOVA
#   write.table(resid(model), paste0("Data/OutputtedData/WSMDP_LinearModel_residuals_",TitleAddendum,"_", carbs[j], ".csv"), col.names=F, row.names=F, sep=",")
# 
#   #calculate model R2
#   r2 <- r.squaredGLMM(model)
#   out <- capture.output(r2[2])
#   cat(paste("The r^2 value of this model for", carbs[j] ,"is", out), file= statsfile, sep="\n", append=TRUE)
#     # Calculate hertiability 
#   model_variances <- as.data.frame(VarCorr(model))
#   #broad sence heritability calculated as variance of genotype/ (variance of geno + var of gxe / num of enviornments + var of error/num of replicates and environ)
#   h2 <- model_variances$vcov[2]/(model_variances$vcov[2]+(model_variances$vcov[1]/4)+(model_variances$vcov[8]/8))
#   out <- capture.output(h2)
#   cat(paste("The heritability of", carbs[j] ,"is", out), file= statsfile, sep="\n", append=TRUE)
# 
#   #Compare model with one with Endosperm
#   cat("Model Comparison, with and without Endosperm Term", file=statsfile, sep="\n", append=TRUE)
#   #####Modeling
#   modelpasteendo<-  paste0(carbs[j], " ~ (1|endo/BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")
#   # Use check if there enough checks in the model set
#   if(length(unique(SampleDFtoModel$Check))>1){
#   modelpasteendo<-  paste0(carbs[j], " ~ Check + (1|endo/BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")
#   }
#   cat(modelpasteendo, file=statsfile, sep="\n", append=TRUE)
#   
#   modelendo <- lmer(modelpasteendo,
#                 data=SampleDFtoModel, REML = TRUE)
#   # Decreasing stopping tolerances
#   strict_tol <- lmerControl(optCtrl=list(xtol_abs=1e-8, ftol_abs=1e-8))
#   if (all(modelendo@optinfo$optimizer=="nloptwrap")) {
#     modelendo <- update(modelendo, control=strict_tol)
#   }
#   # Summary of random effects
#   summary <- summary(modelendo, correlation=FALSE)
#   out <- capture.output(summary)
#   cat(out, file=statsfile, sep="\n", append=TRUE)
#   
#   }
#   
#   # Write out BLUPs for Genotypes
#   write.table(blupHolder, file=paste0("Data/OutputtedData/WSMDP_CarbPheno_InbredBLUPS_",TitleAddendum,".txt"), col.names=T, row.names=F, sep=",")
#   
#   return(blupHolder)
# }
# 
# 
# WFBlups <- linearmodel(CleanedInfoWF,"CleanedOutliersWF")
# NFBlups <- linearmodel(CleanedInfoNF,"CleanedOutliersNF")
# 
# CleanedInfoHWSP_wexcess <- CarbOutlierCleanup(HWSPs,"HWSP",alpha = 0.05)
# CleanedInfoLWSPNF_wexcess <- CarbOutlierCleanup(LWSPNFs,"LWSPNFs",alpha = 0.05)
# CleanedInfoLWSPWF_wexcess <- CarbOutlierCleanup(LWSPWFs,"LWSPWFs",alpha = 0.05)
# 
# CleanedInfoHWSP <- subset(CleanedInfoHWSP_wexcess, !is.na(IsExperimental))
# CleanedInfoLWSPNF <- subset(CleanedInfoLWSPNF_wexcess, !is.na(IsExperimental))
# CleanedInfoLWSPWF <- subset(CleanedInfoLWSPWF_wexcess, !is.na(IsExperimental))
# HWSPBlups <- linearmodel(CleanedInfoHWSP,"CleanedOutliersHWSP")
# LWSPNFBlups <- linearmodel(CleanedInfoLWSPNF,"CleanedOutliersLWSPNF")
# LWSPWFBlups <- linearmodel(CleanedInfoLWSPWF,"CleanedOutliersLWSPWF")



#########################
###Genomic Info###
#########################
# genoinfo <- read.csv("Data/WSMDP_Inbreds.txt",head = FALSE)
genoinfo <- read.csv("Data/WSMDP_Inbreds_2021.6.30.csv",head = TRUE)
# colnames(genoinfo) <- c("Variety","","","","","","","GenoName","source","endo","notes","Region","Program")
colnames(genoinfo) <- c("Index","Variety","Planting2019","Planting20142015","SCMVTest","CAP","Sugar2019","GBS","endo","GenoName")
test1 <- which(genoinfo$Planting20142015==1)
test2 <- which(!is.na(genoinfo$GenoName))
intersect(test1,test2)

#########################
###SNP Relate Establishment###
#########################
infilename <- "WSMDP_SeqG_PreLD"


#read in VCF. This file has been previously filtered for 90% SNP call rate and MinorAlleleFrequency of 0.025
infilename
vcfpath <- paste("Data/RawData/",infilename,".vcf",sep = "")
snpgdsVCF2GDS(vcfpath, "Data/test.gds", method = "biallelic.only")
snpgdsSummary("Data/test.gds")
PreLDGeno <- snpgdsOpen("Data/test.gds")

# #LDVisualization
# L1 <- snpgdsLDMat(PreLDGeno, method="r", slide=250)
# # plot
# LDMatFile <- paste0("Figures/",Sys.Date(),"SNPRelate_LDMatrix.png")
# png(LDMatFile, width = 1000, height = 1000)
# image(abs(L1$LD), col=terrain.colors(64))
# dev.off()

#Lets prune anything with an LD over 0.98
PostLDGeno <- snpgdsLDpruning(PreLDGeno, ld.threshold = 0.98, start.pos = "first", verbose = TRUE)
names(PostLDGeno)
head(PostLDGeno$chr1)
PostLDGeno.id <- unlist(PostLDGeno)

#########################
### PCA ###
#########################
#Read in the SCMV gdsobj with only the snp.id's found by LD pruning and only the sample.id's that we phenotypes

PCA <- snpgdsPCA(PreLDGeno, snp.id = PostLDGeno.id)
#cut off the numbers at the end of the sample.id that don't mean things to humans
holdtrunc<-gsub(PCA$sample.id, pattern = ":.*", replacement = "")
#only include the samples that are in genoinfo
#this is a round about way of looking at only the genotype info I have for the things planted in this trial
#eliminate inbreds tested in other experiments
genoinfo2014only <- genoinfo[which(genoinfo$Planting20142015 == 1),]
#look at what inbreds have the same genoname in both the seq file and the data summary
matching <- data.frame("GenoName" = unique(PCA$sample.id[na.omit(match(genoinfo2014only$GenoName, holdtrunc))]))
#hold the names of the inbreds I have seq info for
holdnocol <- unique(holdtrunc[na.omit(match(genoinfo2014only$GenoName, holdtrunc))])
#combine the geno name with the endosperm mutant time from the genoinfo file
matching$endo <- genoinfo2014only$endo[match(holdnocol,genoinfo2014only$GenoName)]

#conduct pca again with only those samples
PCA <- snpgdsPCA(PreLDGeno, sample.id = matching$GenoName, snp.id = PostLDGeno.id)
pc.perc <- PCA$varprop*100
head(round(pc.perc,2))



#visualize the PCA with ld pruned snps and only the tested samples 
PCAFigureCreation(PCA,pc.perc,matching,infilename,"endo")
# PCAFigureCreation(PCA,pc.perc,matching,filename,"Program")
# PCAFigureCreation(PCA,pc.perc,matching,filename,"Region")


#########################
### Close Snp Relate ###
#########################
outfilename <- "WSMDP_SeqG"
outfile <- paste("Data/RawData/",outfilename ,sep="")
snpgdsGDS2PED(PreLDGeno, outfile, snp.id = PostLDGeno.id)
snpgdsClose(PreLDGeno)




#read in genetic info post MAF 
Seq <- "SeqG"
hmppath <- paste0("Data/RawData/WSMDP_",Seq,".hmp.txt")


# hmppath <- "Data/RawData/WSMDP_SCMV_SeqB.hmp.txt"
# hmppath <- "Data/RawData/WSMDP_SeqC.hmp.txt"
# geno <- fread(hmppath,skip = "rs#")
# str(geno)
# colnames(geno)<-gsub(colnames(geno), pattern = ":.*", replacement = "")
# str(geno)
geno <- read.csv("Data/RawData/SeqG_numericFormat.csv")

BlupGenoCleanup <- function(BlupDFName){
BlupDF <- read.csv(file=paste0("Data/OutputtedData/WSMDP_CarbPheno_InbredBLUPS_CleanedOutliers",BlupDFName,".txt"))#, col.names=T, row.names=F)
  
BlupDFGeno <- merge(BlupDF,genoinfo, by = "Variety")
c1 <- which(colnames(BlupDFGeno)=="Starch.BLUP")
c2 <- which(colnames(BlupDFGeno)=="Total.Sugar.BLUP")
c3 <- which(colnames(BlupDFGeno)=="GenoName")
c4 <- which(colnames(BlupDFGeno)=="endo")

#remove rows with no geno information and remove rows that are duplicates. 
#TODO:Check which duplicates are being deleted and if its a good choice
# BlupDFGeno <- BlupDFGeno[-which(BlupDFGeno$GenoName == ""),]
# BlupDFGeno <- BlupDFGeno[-which(BlupDFGeno$GenoName == "0"),]
BlupDFGeno <- BlupDFGeno[-which(is.na(BlupDFGeno$GenoName)),]
if(length(which(is.na(BlupDFGeno$endo)))>0){
BlupDFGeno <- BlupDFGeno[-which(is.na(BlupDFGeno$endo)),]}
BlupDFGeno <- BlupDFGeno %>% distinct(GenoName, .keep_all = TRUE)

#seperate out superflous information
BlupDFGenoJustPheno <- BlupDFGeno[,c(c3,c1:c2,c4)]
BlupDFGenoJustPheno$endo <- as.factor(BlupDFGenoJustPheno$endo)
return(BlupDFGenoJustPheno)

}



WFBlupsGenoJustPheno <- BlupGenoCleanup("WF")
NFBlupsGenoJustPheno <- BlupGenoCleanup("NF")
HWSPBlupsGenoJustPheno <- BlupGenoCleanup("HWSP")
LWSPNFBlupsGenoJustPheno <- BlupGenoCleanup("LWSPNF")
LWSPWFBlupsGenoJustPheno <- BlupGenoCleanup("LWSPWF")



c1 <- which(colnames(WFBlupsGenoJustPheno)=="Starch.BLUP")
# c1 <- which(colnames(WFBlupsGenoJustPheno)=="Glucose.BLUP")
c2 <- which(colnames(WFBlupsGenoJustPheno)=="Total.Sugar.BLUP")
blups <- colnames(WFBlupsGenoJustPheno[c1:c2])

for(blup in blups){
GWASPolyRunner(WFBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"WFBLUP")
GWASPolyRunner(WFBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"WFBLUP","endo","factor")

GWASPolyRunner(HWSPBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"HWSPBlups")
GWASPolyRunner(HWSPBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"WFBLUP","endo","factor")

GWASPolyRunner(LWSPNFBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPNFBlups")
GWASPolyRunner(LWSPNFBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPNFBlups","endo","factor")

GWASPolyRunner(LWSPWFBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPWFBlups")
GWASPolyRunner(LWSPWFBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPWFBlups","endo","factor")

}


#########################
###GapitKinship Matrix###
#########################
#####Altered from Baseggio's code GWAS_sweet_toco_univ_withcovariate_20180202.R
### Genotypes
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")
path.to.res <- "C:/Users/LHislop/Documents/GitHub/WSMDP_Carb/Data/OutputtedData/GapitOut/"
Seq <- "SeqG"
hmppath <- paste0("Data/RawData/WSMDP_",Seq,".hmp.txt")
geno <- fread(hmppath,skip = "rs#")
str(geno)
colnames(geno)<-gsub(colnames(geno), pattern = ":.*", replacement = "")
colnames(geno)<-gsub(colnames(geno), pattern = "rs#", replacement = "rs")
#assemble the info about marker, chromosome and position

#convert from hapmap to numeric

geno_clean <- geno
geno_clean[which(geno == "A")] <- "AA"
geno_clean[which(geno == "C")] <- "CC"
geno_clean[which(geno == "G")] <- "GG"
geno_clean[which(geno == "T")] <- "TT"
geno_clean[which(geno == "R")] <- "AG"
geno_clean[which(geno == "Y")] <- "CT"
geno_clean[which(geno == "S")] <- "CG"
geno_clean[which(geno == "W")] <- "AT"
geno_clean[which(geno == "K")] <- "GT"
geno_clean[which(geno == "M")] <- "AC"
geno_clean[which(geno == "N")] <- NA
geno_clean[which(geno == "0")] <- NA
geno_clean[which(geno == "-")] <- NA

#some of the SNPs have more than 2 homozygous genotypes. Get rid of them 
toomanyHomos <- which(apply(geno_clean, 1, function(x) length(unique(x)))>4)
if(length(toomanyHomos) > 0){
  geno_clean <- geno_clean[-toomanyHomos,]
  }

str(geno)
nbtraits <- 7
# pheno <- read.table('sweetcorn_toco384_nodent_noHi_NA_20180715.txt', header=T, sep="\t", colClasses=c('factor',rep('numeric', nbtraits)), as.is=T)
nlines <- nrow(WFBlupsGenoJustPheno)

# geno <- read.table(hmppath, header=FALSE, sep = "\t")
# str(geno)
# colnames(geno)<-gsub(colnames(geno), pattern = ":.*", replacement = "")
str(geno)
setwd(path.to.res)

kin.alg <- "VanRaden"
res <- GAPIT( G = geno, kinship.cluster = c("average"), kinship.group = c("Mean"), kinship.algorithm = kin.alg,SNP.impute = "Middle", Major.allele.zero = TRUE)


############################################
####### Run GWAS for TOCOCHROMANOLS ########
############################################

## Read the files
# setwd("/local/workdir/mb2446/GWAS/")


covariates <- rbind.data.frame(c('Starch','Starch'),
                               c('Total.Polysaccharides','Total.Polysaccharides'),
                               c('WSP','WSP'),
                               c('Glucose','Glucose'),
                               c('Fructose','Fructose'),
                               c('Sucrose','Sucrose'),
                               c('Total.Sugar','Total.Sugar'))
colnames(covariates) <- c('trait','snp')

# geno <- read.delim('sweet_parimp_toco384_174K.hmp.txt', header=FALSE, sep='\t')

nbtraits <- 7
# pheno <- read.table('sweetcorn_toco384_nodent_noHi_NA_20180715.txt', header=T, sep="\t", colClasses=c('factor',rep('numeric', nbtraits)), as.is=T)
# path.to.res <- "/local/workdir/mb2446/GWAS/Univariate/Toco/0PC.2yr.K11.174K.384taxa.middle.na.KernelType.20180715/"

# parent.dir <- "/local/workdir/mb2446/GWAS/"
pheno <- WFBlupsGenoJustPheno

# library(beepr)
# system.time(beep() sound = "fanfare"))
  for (i in seq(1,7)) {
    #### Current trait ####
    curr.trait <- colnames(pheno[ (1 + i) ])
    # curr.trait <- as.character(covariates[i,1])
    # system(paste('mkdir -p ', path.to.res, curr.trait, sep=''))
    setwd(paste0(path.to.res ,curr.trait))
    j <- which(colnames(pheno)==curr.trait)
    curr.pheno = pheno[,c(1,j)]
    #### User input Kinship
    myKI <- read.table(paste(path.to.res, '/GAPIT.Kin.VanRaden.csv', sep=''), header=FALSE, sep=',')
    # myCV <- read.table(paste('/local/workdir/mb2446/GWAS/Covariates/KernelType.txt', sep=''), header=T, sep='\t')
    res <- GAPIT( Y = curr.pheno, G = geno, KI=myKI, PCA.total=3, plot.style="Ocean")
  }
