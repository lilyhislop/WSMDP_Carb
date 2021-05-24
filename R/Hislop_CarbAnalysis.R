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
library(GWASpoly)
library("bigmemory") #to make a matrix big
# library(rrBLUP)

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
source("R/hmpToNumeric.R")


#Read in the sample information.
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
  return(OutDF)
}

HWSPsDF <- InfoCombination(HWSPsDF,SampleInfo,BookInfo)
LWSPWFsDF <- InfoCombination(LWSPWFsDF,SampleInfo,BookInfo)
LWSPNFDF <- InfoCombination(LWSPNFDF,SampleInfo,BookInfo)


#Eliminate irrelivant endosperm mutants from each df. We don't want to include sh2 samples that were predicted by the HWSP calibrated equations
#should sh2i samples be estimated by the sh2 calibrated equations or the other calibrated equation?
HWSPs <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se" | (is.na(HWSPsDF$endo)& (HWSPsDF$Year== "15"|HWSPsDF$Year == "14"))),]
# HWSPs <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se" | is.na(HWSPsDF$endo)),]
LWSPWFs <- LWSPWFsDF[which(LWSPNFDF$endo != "se" & LWSPNFDF$endo != "su1"),]
LWSPNFs <- LWSPNFDF[which(LWSPNFDF$endo == "sh2" | LWSPNFDF$endo == "sh2i"),]

#lets look at the validation data that I predicted from kahtleen and Jareds data and output it for further analysis
jared <- HWSPsDF[which(HWSPsDF$Year == "sc"),]
write.csv(jared, file = "Data/OutputtedData/JaredsPredictedLineswHWSPeqn.csv")
kathleen <- HWSPsDF[which(HWSPsDF$Year != "13"&HWSPsDF$Year != "14"&HWSPsDF$Year != "15"&HWSPsDF$Year != "sc"),]
write.csv(kathleen,file = "Data/OutputtedData/KathleenPredictedLineswHWSPeqn.csv")

########mash it all together ######
#All the predictions with the field lines
CarbInfoExpandedWFDF <- rbind(HWSPs,LWSPWFs)
#all the predictions without the field lines
CarbInfoExpandedNFDF <- rbind(HWSPs,LWSPNFs)

#add in a year/location environment factor
CarbInfoExpandedWFDF$Envi <- paste(CarbInfoExpandedWFDF$Year,CarbInfoExpandedWFDF$Location)
CarbInfoExpandedNFDF$Envi <- paste(CarbInfoExpandedNFDF$Year,CarbInfoExpandedNFDF$Location)

#visualize these data sets pre cleaning
CarbDataFrameVis(CarbInfoExpandedWFDF,"WithField_WithOutliers")
CarbDataFrameVis(CarbInfoExpandedNFDF,"NoField_WithOutliers")

#clean up the predictive data frames. Reassign or delete outliers 
CleanedInfoWF <- CarbOutlierCleanup(CarbInfoExpandedWFDF,"WF",alpha = 0.05)
CleanedInfoNF <- CarbOutlierCleanup(CarbInfoExpandedNFDF,"NF",alpha = 0.05)
write.csv(file = "Data/OutputtedData/CleanedInfoWFOutput.csv",CleanedInfoWF)
write.csv(file = "Data/OutputtedData/CleanedInfoNFOutput.csv",CleanedInfoNF)
LotsoLabels <- 

#revisualize the dataframes 
CarbDataFrameVis(CleanedInfoWF,"WithField_Cleaned")
CarbDataFrameVis(CleanedInfoNF,"NoField_Cleaned")

#write the names of the varieties used to a csv file so we can find the corresponding GBS data
write.csv(file = "Data/OutputtedData/InbredsWithinWSMDPCarbDataNF.csv",unique(CleanedInfoNF[c("Variety", "endo")]))
write.csv(file = "Data/OutputtedData/InbredsWithinWSMDPCarbDataWF.csv",unique(CleanedInfoWF[c("Variety", "endo")]))


#######Equation Validation!##############
#Now I have a variable that has all the projected values, for the values used to calibrate the equations. WHat are the statistics on that?
#set carb to a number 1:7. carb <- c(Fructose, Glucose, Sucrose, Total Sugar, Starch, Total Polysaccharide, WSP)
#this functions needs to have the data frame set up so the carbs being compared are directly next to eachother 
EqnStats <- function(DF){
  
  #establish dataframe used to record stats
  Out <- data.frame(Carb = c("Fructose", "Glucose", "Sucrose", "Total Sugar", "Starch", "Total Polysaccharide", "WSP"),
                    RMSEP = rep(NA,7), 
                    bias = rep(NA,7),
                    SEE = rep(NA,7),
                    slope = rep(NA,7),
                    intercept = rep(NA,7),
                    R2 = rep(NA,7))
  dfpos <- c(2,4,6,8,10,12,14)
  for(carb in 1:7){
    #Calculated the RMSEP
    Out$RMSEP[carb] <- sqrt(sum((DF[,dfpos[carb]]- DF[,dfpos[carb]-1])^2, na.rm = TRUE)/dim(DF)[1])
    #Calculate the Bias
    Out$bias[carb] <- mean(DF[,dfpos[carb]], na.rm = TRUE) - mean(DF[,dfpos[carb]-1],na.rm = TRUE)
    #Calculated the SEE
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


######Visualize the validation. Look at all the samples that were wetlabbed and used to calibrate the equation
wetlab <- read.csv("Data/WSMDP_Wetlab_StarchSugarData_FormatedForWinISI_WithR.csv")
wetlabDF <- wetlab[,c(1,6:8,3,2,4,5)]

WLLabels <- c("Samples","Starch_WL","Total.Polysaccharides_WL", "WSP_WL","Glucose_WL","Fructose_WL","Sucrose_WL","Total.Sugar_WL")
colnames(wetlabDF) <- WLLabels 

WFWL <- merge(CleanedInfoWF, wetlabDF, by = "Samples")
NFWL <- merge(CleanedInfoNF, wetlabDF, by = "Samples")
alternatingorder <- c(1,5,28,6,29,7,30,8,31,9,32,10,33,11,34)
WFWLdf <- WFWL[,alternatingorder]
NFWLdf <- NFWL[,alternatingorder]
WFWLdfEqnStatsR <- R2Vis(WFWLdf[,2:15], "CleanedWSPeqnWF_PredVsWetlab_for_Calibration_Samples", EqnStats(WFWLdf[,2:15]))
NFWLdfEqnStatsR <- R2Vis(NFWLdf[,2:15], "CleanedWSPeqnNF_PredVsWetlab_for_Calibration_Samples", EqnStats(NFWLdf[,2:15]))

write.csv(WFWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsAllWetlabWF.csv")
write.csv(NFWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsAllWetlabNF.csv")

######Visualize the validation. Only the samples that were NOT used to create the equation
valwetlab <- read.csv("Data/WSMDP_EqnValidation_Wetlab_Data.csv")
valwetlabDF <- valwetlab[-2]
colnames(valwetlabDF)[1] <- "Samples"

WFValWL <- merge(CleanedInfoWF, valwetlabDF, by = "Samples")
NFValWL <- merge(CleanedInfoNF, valwetlabDF, by = "Samples")
WFValWLdf <- WFValWL[,alternatingorder]
NFValWLdf <- NFValWL[,alternatingorder]
WFValWLdfEqnStatsR <- R2Vis(WFValWLdf[,2:15], "CleanedWSPeqnWF_PredVsWetlab_ValidationSubset", EqnStats(WFValWLdf[,2:15]))
NFValWLdfEqnStatsR <- R2Vis(NFValWLdf[,2:15], "CleanedWSPeqnNF_PredVsWetlab_ValidationSubset", EqnStats(NFValWLdf[,2:15]))

write.csv(WFValWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationWF.csv")
write.csv(NFValWLdfEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationNF.csv")


####With Jared WetlabDATA
HWSPsJ <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se" | is.na(HWSPsDF$endo)),]
alternatingorder2 <- c(1,5,26,6,27,7,28,8,29,9,30,10,31,11,32)
CarbInfoExpandedWFJDF <- rbind(HWSPsJ,LWSPWFs)
CarbInfoExpandedNFJDF <- rbind(HWSPsJ,LWSPNFs)
WFValWLJ <- merge(CarbInfoExpandedWFJDF, valwetlabDF, by = "Samples")
NFValWLJ <- merge(CarbInfoExpandedNFJDF, valwetlabDF, by = "Samples")
WFValWLJdf <- WFValWLJ[,alternatingorder2]
NFValWLJdf <- NFValWLJ[,alternatingorder2]
WFValWLdfJEqnStatsR <- R2Vis(WFValWLJdf[,2:15], "UnCleanedWSPeqnWF_PredVsWetlab_ValidationSubset_WJared", EqnStats(WFValWLdf[,2:15]))
NFValWLdfJEqnStatsR <- R2Vis(NFValWLJdf[,2:15], "UnCleanedWSPeqnNF_PredVsWetlab_ValidationSubset_WJared", EqnStats(NFValWLdf[,2:15]))

write.csv(WFValWLdfJEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationWF_WJared.csv")
write.csv(NFValWLdfJEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationNF_WJared.csv")



##########Linear Model Analysis!##########
head(CleanedInfoWF)

# Which model is best? function
ModelCheck <- function(DF){
  bestAIC <- 10000000
  bestModel <- "Test"
# first, lets clean up the data so its readable and good
DFSubset <- subset(DF, !is.na(IsExperimental))
str(DFSubset)
colnames(DFSubset)

DFSubset$superblock <- as.factor(DFSubset$superblock)
DFSubset$Col <- as.factor(DFSubset$Col)
DFSubset$Row <- as.factor(DFSubset$Row)
DFSubset$Year <- as.factor(DFSubset$Year)
DFSubset$Check <- as.factor(DFSubset$Check)
DFSubset$block <- as.factor(DFSubset$block)
DFSubset$Rep <- as.factor(DFSubset$Rep)
DFSubset$PlotNum <- as.factor(DFSubset$PlotNum)

# Firstmodel: Most basic

# ModelA
# G, E, GxE, row and column in environment
#with block, superblock as stand alone fixed effects, envi as random
modelA <- lmer(Total.Sugar ~ Check + (1|Variety) + Envi + (1|Envi/Row) + (1|Envi/Col),data=DFSubset, REML = TRUE)
AIC(modelA)
if(AIC(modelA) < bestAIC){bestAIC = AIC(modelA)
bestModel = "A"}

# ModelB
#with block, superblock, row, column as stand alone fixed effects, envi as random
modelB <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + Row + Col + block, data=DFSubset, REML = TRUE)
AIC(modelB)
if(AIC(modelB) < bestAIC){bestAIC = AIC(modelB)
bestModel = "B"}

# ModelC
#with block, superblock as stand alone fixed effects, envi as random
modelC <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + block, data=DFSubset, REML = TRUE)
AIC(modelC)
if(AIC(modelC) < bestAIC){bestAIC = AIC(modelC)
bestModel = "C"}

# Model D
#with block, superblock, row, column as stand alone fixed effects, envi as random
modelD <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + Row + Col + block, data=DFSubset, REML = TRUE)
AIC(modelD)
if(AIC(modelD) < bestAIC){bestAIC = AIC(modelD)
bestModel = "D"}

# Model E
#with block, superblock, row,  as stand alone fixed effects, envi as random
modelE <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + Row +  block, data=DFSubset, REML = TRUE)
AIC(modelE)
if(AIC(modelE) < bestAIC){bestAIC = AIC(modelE)
bestModel = "E"}

# Model F
#with block, superblock, row, column  effects and envi all as random
modelF <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col) , data=DFSubset, REML = TRUE)
AIC(modelF)
if(AIC(modelF) < bestAIC){bestAIC = AIC(modelF)
bestModel = "F"}

# Model G
#with block, superblock nested random effect
modelG <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi/superblock/block), data=DFSubset, REML = TRUE)
AIC(modelG)
if(AIC(modelG) < bestAIC){bestAIC = AIC(modelG)
bestModel = "G"}

# Model H
#model worked out with raegan
modelH<- lmer(Total.Sugar ~  (1|superblock:Variety) + (1+Envi|Variety) + (1|superblock/block), 
                data=DFSubset, REML = TRUE)
AIC(modelH)
if(AIC(modelH) < bestAIC){bestAIC = AIC(modelH)
bestModel = "H"}

# Model J
#model worked out with raegan
modelJ<- lmer(Total.Sugar ~  Check + (1|Variety) + (1|Variety:Envi) + superblock + (1|Envi/superblock/block), 
              data=DFSubset, REML = TRUE)
AIC(modelJ)
if(AIC(modelJ) < bestAIC){bestAIC = AIC(modelJ)
bestModel = "J"}

# Model I
#model worked out with raegan
modelI<- lmer(Total.Sugar ~  Check + (1|Envi/Variety) + (1|superblock) + (1|superblock:Envi) + (1|Envi/block), 
              data=DFSubset, REML = TRUE)
AIC(modelI)
if(AIC(modelI) < bestAIC){bestAIC = AIC(modelI)
bestModel = "I"}

# Model K. This is the one from the stats adviser
modelK<- lmer(Total.Sugar ~  Check + Envi + (1|Envi:Variety) + (1|superblock/Variety) + superblock + (1|superblock/block), 
              data=DFSubset, REML = TRUE)
AIC(modelK)
if(AIC(modelK) < bestAIC){bestAIC = AIC(modelK)
bestModel = "K"}

return(bestModel)

}



ModelCheck(CleanedInfoWF)
ModelCheck(CleanedInfoNF)

# best model is B
modelB <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + Row + Col + block, data=CleanedInfoWF, REML = TRUE)
AIC(modelB)
summary(modelB)
plot(modelB)
qqmath(ranef(modelB, condVar=TRUE))
#looks pretty ok, probably don't need tranformation
qqnorm(residuals(modelB))


# Compare this model to the one the Baseggio put in his paper
modelBaseggio <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + (1|Envi/superblock) + (1|Envi/Row) + (1|Envi/Col) + (1|Envi/superblock/block), data=CleanedInfoWF, REML = TRUE)
AIC(modelBaseggio)
summary(modelBaseggio)
plot(modelBaseggio)
qqmath(ranef(modelBaseggio, condVar=TRUE))
#looks pretty ok, probably don't need tranformation
qqnorm(residuals(modelBaseggio))

anova(modelBaseggio, modelB)

#ok, lets go with matts model


#establish dataframe to store variances from each factor
VarDF <- data.frame("Carb" = colnames(CleanedInfoWF)[5:11],"Variety" = rep(NA,7),"Envi" = rep(NA,7),"Variety:Envi" = rep(NA,7),  "superblock" = rep(NA,7), "block" = rep(NA,7),"Col" = rep(NA,7),"Row" = rep(NA,7),"Residuals"= rep(NA,7))

#for each carb version, look what factors are incluencing the variation
for(i in 1:7){
  #formula is carb ~ Gene + Envi + Gene*Envi + rep + endosperm type + error. Should add block in as well 
  formula1 <- paste(colnames(CleanedInfoWF)[i+4],"~ Variety*Envi + Rep + endo")
  fit1 <- lm(formula1,data=CleanedInfoWF)
AIC1 <- extractAIC(modelBaseggio)
AIC1
summary(modelBaseggio)$r.square
print(anova(modelBaseggio))
out <- anova(modelBaseggio)
SStotal <- sum(out$`Sum Sq`)
for(j in 1:6){
  #variance explained is caluclated by the sum of squares divided by the sum of squares total
VarDF[i,j+1] <- (out$`Sum Sq`[j]/SStotal)}
}
VarDFMelt <- melt(VarDF)



#######Graph the different variances explained by different factors######
png(paste("Figures/WSMDP_AllNIRPred_MixedEqn_PercentVarianceExplainedby_Factors_WithField.png",sep=""), width = 1000, height = 500)
barchart(~value|variable, group = factor(Carb), data= VarDFMelt,main = "Percent Phenotypic Variance Explained",layout = c(6,1),
         key = simpleKey(text = colnames(CleanedInfoWF)[5:11],
                         rectangles = TRUE, points = FALSE, space = "right"))
dev.off()

#########Some Other Stats########
carbs <- colnames(CleanedInfoWF)[5:11]
h2DF <- data.frame("Carb" = colnames(CleanedInfoWF)[5:11],"mean" = c(1:7), "SD" = c(1:7),"min" = c(1:7),"max" = c(1:7),"h2" = c(1:7))
#establish dataframe to store variances from each factor
VarDF <- data.frame("Carb" = colnames(CleanedInfoWF)[5:11],"Variety" = rep(NA,7),"Envi" = rep(NA,7), "Rep" = rep(NA,7),"endo" = rep(NA,7),"Variety:Envi" = rep(NA,7),"Residuals"= rep(NA,7))



for(carb in carbs){
normality <- shapiro.test(pull(CleanedInfoWF[,carb]))
  out <- capture.output(normality)


# model <- lmer(get(carb) ~ (1|Variety) + (1|Envi) + (1|Envi/Rep)+ (1|Variety:Envi),data=CleanedInfoWF, REML = TRUE)
summary(modelBaseggio, correlation = FALSE)
random_effects <- ranef(modelBaseggio)

#Write out the blups
write.table(random_effects$Variety, paste0("Data/OutputtedData/blups_", carb, ".csv"), col.names=F, row.names=F, sep=",")
summary <- summary(modelBaseggio, correlation = FALSE)

#Write out the residuals
write.table(resid(modelBaseggio), paste0("Data/OutputtedData/resids_", carb, ".csv"), col.names=F, row.names=F, sep=",")

#####Code for heritibility taken from jonathan Renk
# Calculate hertiability 
model_variances <- as.data.frame(VarCorr(modelBaseggio))
h2 <- model_variances$vcov[2]/(model_variances$vcov[2]+(model_variances$vcov[1]/4)+(model_variances$vcov[6]/8))
####Record Data
h2DF[which(h2DF$Carb == carb),"h2"] <- h2
h2DF[which(h2DF$Carb == carb),"mean"] <- mean(pull(CleanedInfoWF[,carb]),na.rm = TRUE)
h2DF[which(h2DF$Carb == carb),"min"] <- min(pull(CleanedInfoWF[,carb]),na.rm = TRUE)
h2DF[which(h2DF$Carb == carb),"max"] <- max(pull(CleanedInfoWF[,carb]),na.rm = TRUE)
h2DF[which(h2DF$Carb == carb),"SD"] <- SD(pull(CleanedInfoWF[,carb]),na.rm = TRUE)
out <- capture.output(h2)
cat(out, file=paste0("Data/OutputtedData/h2_",carb,".txt"), sep="\n", append=TRUE)

pdf(paste0("Figures/assumptions_", carb, ".pdf"), width = 15, height = 5)
par(mfrow=c(1,2))


# Model Fit with REML
plot(fitted(modelBaseggio), residuals(modelBaseggio), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
abline(h=0,col="red", lwd=1, lty=1)
# histogram of residuals
hist(residuals(modelBaseggio),main="Histogram of residuals",freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
x=seq(-5e-15,9e-15,5e-15)
curve(dnorm(x,mean(residuals(modelBaseggio)),sd(residuals(modelBaseggio))),add=T,lwd=2, col="red", lty=1)
# # qq plot
# qqPlot(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles") 

dev.off()


}

write.csv(h2DF, file=paste0("Data/OutputtedData/h2_ALLCarb.txt"), append=FALSE)

#Number of total Samples
str(unique(CleanedInfoWF$Samples))
#Number of varieties included
str(unique(CleanedInfoWF$Variety))

########Plot Correlations########
justthebitsNF <- CleanedInfoNF[5:11]
png(paste("Figures/WSMDP_AllNIRPred_MixedEqn_CorrelationfromPSYCH_NoField.png",sep=""), width = 500, height = 500)
pairs.panels(justthebitsNF, scale = TRUE)
dev.off()

justthebitsWF <- CleanedInfoWF[5:11]
png(paste("Figures/WSMDP_AllNIRPred_MixedEqn_CorrelationfromPSYCH_WithField.png",sep=""), width = 500, height = 500)
pairs.panels(justthebitsWF, scale = TRUE)
dev.off()


########GWAS ZONE##########
##First things first in the GWAS zone. Find the GBS Names associated with the variety
# vcfFilename <- "Data/RawData/WSMDP_SCMV_SeqE_Vcf.vcf"
# vcf <- readVcf(vcfFilename)
# snpmat <- genotypeToSnpMatrix(vcf)
# Matrix <- snpmat$genotypes@.Data



genoinfo <- read.csv("Data/WSMDP_Inbreds.txt",head = FALSE)
colnames(genoinfo) <- c("Variety","","","","","","","GenoName","source","endo","notes","Region","Program")
myG <- read.delim("Data/RawData/WSMDP_SCMV_SeqE.hmp.txt", head = FALSE)
myGnames <- myG[1,]
myGnames<-gsub(myGnames, pattern = ":.*", replacement = "")
myG[1,] <- myGnames



# 
# CleanedInfoNFwGeno <- merge(CleanedInfoNF,genoinfo, by = "Variety")
# myY <- CleanedInfoNFwGeno[,c(19,6:12)]#,"class","tbl")
# ######GWAS ALL TOGETHER #######
# setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb/Data/OutputtedData/GAPIT")
# #A plain GAPIT MLM
# # system.time(
# # myGAPIT <- GAPIT(
# #     Y=myY,
# #     G=myG,
# #     PCA.total=3
# #     ))
#   # FarmCPU
# system.time(
# myGAPIT <- GAPIT(
#   Y=myY,
#   G=myG,
#   PCA.total=3,
#   model = "FarmCPU"
# )
# )
# # MLMM & gBlup
# system.time(
#   myGAPIT <- GAPIT(
#     Y=myY,
#     G=myG,
#     PCA.total=3,
#     model = c("MLMM","gBLUP")
#   )
# )
#GAPIT to analyze phenotype info
# myPhenotypes <- GAPIT.Phenotype.View(myY = myY,)

######Function to look at HWSP, LWSP subsets separateluy######
GAPITRunner <- function(DF, label,Seq){
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")
CleanedDF <- CarbOutlierCleanup(DF,label,alpha = 0.05)
CleanedDFwGeno <- merge(CleanedDF,genoinfo, by = "Variety")
c1 <- which(colnames(CleanedDFwGeno)=="Starch")
c2 <- which(colnames(CleanedDFwGeno)=="Total.Sugar")
c3 <- which(colnames(CleanedDFwGeno)=="GenoName")
myY <- CleanedDFwGeno[,c(c3,c1:c2)]
myY <- myY[-which(myY$GenoName == ""),]
myYMerged <- myY %>%
  group_by(GenoName) %>%
  summarise(Starch = mean(Starch,na.rm = T),Total.Polysaccharides = mean(Total.Polysaccharides,na.rm = T), 
            WSP = mean(WSP,na.rm = T), Glucose =mean(Glucose,na.rm = T), Fructose = mean(Fructose,na.rm = T),
            Sucrose = mean(Sucrose,na.rm = T), Total.Sugar = mean(Total.Sugar,na.rm = T))
  
setwd(paste("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb/Data/OutputtedData/GAPIT/",label,"/",Seq,sep = ""))
myGAPIT <- GAPIT(
    Y=myY,
    G=myGB,
    PCA.total=3,
    model = c("MLMM","gBLUP"))#,"FarmCPU")
# myGAPIT <- GAPIT(
#     Y=myY,
#     G=myG,
#     PCA.total=5,
#     method.bin = "optimum",
#     model = "FarmCPU"
#   
# )
}
system.time({
GAPITRunner(CarbInfoExpandedWFDF, "WF","SeqB")
GAPITRunner(CarbInfoExpandedNFDF, "NF","SeqB")
GAPITRunner(HWSPs,"HWSP","SeqB")
GAPITRunner(LWSPWFs,"LWSPWF","SeqB")
GAPITRunner(LWSPNFs,"LWSPNF","SeqB")})


CleanedDFwGeno <- merge(CleanedInfoWF,genoinfo, by = "Variety")
c1 <- which(colnames(CleanedDFwGeno)=="Starch")
c2 <- which(colnames(CleanedDFwGeno)=="Total.Sugar")
c3 <- which(colnames(CleanedDFwGeno)=="GenoName")
myY <- CleanedDFwGeno[,c(c3,c1:c2)]
myY <- myY[-which(myY$GenoName == ""),]
myYMerged <- myY %>%
  group_by(GenoName) %>%
  summarise(Starch = mean(Starch,na.rm = T),Total.Polysaccharides = mean(Total.Polysaccharides,na.rm = T), 
            WSP = mean(WSP,na.rm = T), Glucose =mean(Glucose,na.rm = T), Fructose = mean(Fructose,na.rm = T),
            Sucrose = mean(Sucrose,na.rm = T), Total.Sugar = mean(Total.Sugar,na.rm = T))
post <- data.frame(myYMerged)


setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb/Data/OutputtedData/GAPIT/WF/FarmCPU")
myGAPIT <- GAPIT(
  G = myG, output.numerical = TRUE
)
myGD <- read.big.matrix("GAPIT.Genotype.Numerical.txt", type="char", sep="\t", head = TRUE)
myGM <- read.table("GAPIT.Genotype.map.txt", head = TRUE)

myGAPIT <- GAPIT( 
  Y=post[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  method.bin="optimum",
  model="FarmCPU"
)

pvals <- FarmCPU.P.Threshold(
  Y=myYMerged[,c(1,2)], #only two columns allowed, the first column is taxa name and the second is phenotype
  GD = myGD,
  GM = myGM,
  trait="Starch", #name of the trait, only used for the output file name
  theRep=100 #number of permutation times 
)


######with GenoSeq B#########
myGB <- read.delim("Data/RawData/WSMDP_SCMV_SeqB.hmp.txt", head = FALSE)
myGBnames <- myGB[1,]
myGBnames<-gsub(myGBnames, pattern = ":.*", replacement = "")
myGB[1,] <- myGBnames


CleanedDFwGeno <- merge(CleanedInfoWF,genoinfo, by = "Variety")
c1 <- which(colnames(CleanedDFwGeno)=="Starch")
c2 <- which(colnames(CleanedDFwGeno)=="Total.Sugar")
c3 <- which(colnames(CleanedDFwGeno)=="GenoName")
myY <- CleanedDFwGeno[,c(c3,c1:c2)]
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb/Data/OutputtedData/GAPIT/WF/FarmCPU")
myGAPIT <- GAPIT(
  G = myGB, output.numerical = TRUE
)
myGD <- read.big.matrix("GAPIT.Genotype.Numerical.txt", type="char", sep="\t", head = TRUE)
myGM <- read.table("GAPIT.Genotype.map.txt", head = TRUE)


myFarmCPU<-FarmCPU(Y=myYMerged[,c(1,3)],
                         GD=myGD,
                         GM=myGM)
myGAPIT <- GAPIT( 
  Y=myYMerged[,c(1,2)], 
  GD=myGD,
  GM=myGM,
  PCA.total=5,
  method.bin="optimum",
  model="FarmCPU"
)

pvals <- FarmCPU.P.Threshold(
  Y=myYMerged[,c(1,2)], #only two columns allowed, the first column is taxa name and the second is phenotype
  GD = myGD,
  GM = myGM,
  trait="Starch", #name of the trait, only used for the output file name
  theRep=100 #number of permutation times 
)







#########GWASPOLY RUNNING#########
#read in genetic info post MAF and LD pruning. Pruning done by SNPrelate package. Outputted previously as a plink format, converted to Hapmap by tassel and read back in
hmppath <- "Data/RawData/WSMDP_SCMV_SeqB.hmp.txt"
SCMV_geno <- fread(hmppath,skip = "rs#")
geno_scmv <- SCMV_geno

str(geno_scmv)
colnames(geno_scmv)<-gsub(colnames(geno_scmv), pattern = ":.*", replacement = "")
str(geno_scmv)


#########################
### GWASPoly Function ###
#########################
#trait to analyze
# CleanedDF <- CarbOutlierCleanup(DF,label,alpha = 0.05)
CleanedDFwGeno <- merge(CleanedInfoWF,genoinfo, by = "Variety")
c1 <- which(colnames(CleanedDFwGeno)=="Starch")
c2 <- which(colnames(CleanedDFwGeno)=="Total.Sugar")
c3 <- which(colnames(CleanedDFwGeno)=="GenoName")
c4 <- which(colnames(CleanedDFwGeno)=="endo.x")
myYWF <- CleanedDFwGeno[,c(c3,c1:c2,c4)]
myYWF <- myYWF[-which(myYWF$GenoName == ""),]
myYWFmerged <- myYWF %>%
  group_by(GenoName) %>%
  summarise(Starch = mean(Starch,na.rm = T),Total.Polysaccharides = mean(Total.Polysaccharides,na.rm = T), 
            WSP = mean(WSP,na.rm = T), Glucose =mean(Glucose,na.rm = T), Fructose = mean(Fructose,na.rm = T),
            Sucrose = mean(Sucrose,na.rm = T), Total.Sugar = mean(Total.Sugar,na.rm = T), endo.x = first(endo.x))

CleanedDFwGeno <- merge(CleanedInfoNF,genoinfo, by = "Variety")
c1 <- which(colnames(CleanedDFwGeno)=="Starch")
c2 <- which(colnames(CleanedDFwGeno)=="Total.Sugar")
c3 <- which(colnames(CleanedDFwGeno)=="GenoName")
c4 <- which(colnames(CleanedDFwGeno)=="endo.x")
myYNF <- CleanedDFwGeno[,c(c3,c1:c2,c4)]
myYNF <- myYNF[-which(myYNF$GenoName == ""),]

myYNFmerged <- myYNF %>%
  group_by(GenoName) %>%
  summarise(Starch = mean(Starch,na.rm = T),Total.Polysaccharides = mean(Total.Polysaccharides,na.rm = T), 
            WSP = mean(WSP,na.rm = T), Glucose =mean(Glucose,na.rm = T), Fructose = mean(Fructose,na.rm = T),
            Sucrose = mean(Sucrose,na.rm = T), Total.Sugar = mean(Total.Sugar,na.rm = T), endo.x = first(endo.x))


CleanedDF <- CarbOutlierCleanup(HWSPs,"HWSPs",alpha = 0.05)
CleanedDFwGeno <- merge(CleanedDF,genoinfo, by = "Variety")
c1 <- which(colnames(CleanedDFwGeno)=="Starch")
c2 <- which(colnames(CleanedDFwGeno)=="Total.Sugar")
c3 <- which(colnames(CleanedDFwGeno)=="GenoName")
c4 <- which(colnames(CleanedDFwGeno)=="endo.x")
myYHWSP <- CleanedDFwGeno[,c(c3,c1:c2,c4)]
myYHWSP <- myYHWSP[-which(myYHWSP$GenoName == ""),]
myYHWSPMerged <- myYHWSP %>%
  group_by(GenoName) %>%
  summarise(Starch = mean(Starch,na.rm = T),Total.Polysaccharides = mean(Total.Polysaccharides,na.rm = T), 
            WSP = mean(WSP,na.rm = T), Glucose =mean(Glucose,na.rm = T), Fructose = mean(Fructose,na.rm = T),
            Sucrose = mean(Sucrose,na.rm = T), Total.Sugar = mean(Total.Sugar,na.rm = T), endo.x = first(endo.x))

system.time(
for(carb in carbs){
GWASPolyRunner(myYWFmerged[,1:8],geno_scmv,carb,"NoFixedEffect_PermThresh","SeqB","WFMeaned")
GWASPolyRunner(myYWFmerged,geno_scmv,carb,"EndoFixedEffect_PermThresh","SeqB","WFMeaned","endo.x","factor")
GWASPolyRunner(myYNFmerged[,1:8],geno_scmv,carb,"NoFixedEffect_PermThresh","SeqB","NFMeaned")
GWASPolyRunner(myYNFmerged,geno_scmv,carb,"EndoFixedEffect_PermThresh","SeqB","NFMeaned","endo.x","factor")
GWASPolyRunner(myYHWSPMerged[,1:8],geno_scmv,carb,"NoFixedEffect_PermThresh","SeqB","HWSPMeaned")
GWASPolyRunner(myYHWSPMerged,geno_scmv,carb,"EndoFixedEffect_PermThresh","SeqB","HWSPMeaned","endo.x","factor")

}
)

