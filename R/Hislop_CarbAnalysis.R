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
library(multcomp)
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
# source("http://zzlab.net/GAPIT/gapit_functions.txt")
# source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")
library("devtools")#for intsalling from github
install_github("jendelman/GWASpoly")
library(GWASpoly)#for running Gwas
# library("bigmemory") #to make a matrix big
# library(rrBLUP)
library(emmeans)
library(bigmemory)
library(MuMIn)
library(SNPRelate)
library(beepr)#beep when code is done
library("cAIC4")
library(ggpubr)
library("PAST")
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
source("R/getmedianLDVis.R")
source("R/getpercentileLDVis.R")
# source("R/read.GWASpoly.R")

#########################
###Read in Genomic Info###
#########################
genoinfo <- read.csv("Data/RawData/WSMDP_Inbreds_2021.9.8.csv",head = TRUE)


#########################
###Read in sample data###
#########################

#Read in the sample information.
#this file contains information about the superblock, row, column data for each sample
BookInfo <- read.csv("Data/20142015_WSMDP_PlantingInfo_Condensation.csv")
#There are some duplicates in planting info condensation. Remove duplicates
BookInfo <- BookInfo[!duplicated(BookInfo$SampleLabel),]
length(unique(BookInfo$Inbred))



#this file containes the year, location,Inbred and endosperm of every sample
SampleInfo <- read.csv("Data/RawData/WSMDP_2014-2015_WINY_SampleInfo.csv")
#Fix the sample numbers that are labeled in the NY style e.g " 14A0255" or "15A0005" and change to just the row numbers
SampleInfo$Row <- SampleInfo$Sample.Number
SampleInfo$Row[which(nchar(SampleInfo$Row)>6)] <- substr(SampleInfo$Row[which(nchar(SampleInfo$Row)>6)],4,15)
tail(SampleInfo)
#Generate NIR Code. NIR codes are all YearLocationRow - Rep
SampleInfo$NIRBase <- paste(substr(SampleInfo$year,3,4), ifelse(SampleInfo$location == "Wisconsin", "W","NY"), SampleInfo$Row,sep = "")

#look at the sample info file to make sure it looks ok
tail(SampleInfo)
str(unique(SampleInfo$Variety))
summary.factor(SampleInfo$endo)
length(unique(SampleInfo$Variety))
setdiff(unique(BookInfo$Inbred),unique(SampleInfo$Variety))
sampledInbreds<- intersect(unique(BookInfo$Inbred),unique(SampleInfo$Variety))
SampleInfo$endo <- as.factor(SampleInfo$endo)

#########################
###Read in Wetlab Starch Data###
#########################
###################
wetLabStarch <- read.csv("Data/RawData/WSMDP_Wetlab_Starchs.csv")

#separate out the data that aren't controls
wetLabStarchNoCtrl<- wetLabStarch[which(wetLabStarch$Sample_ID != "control"),c(5,9:11)]
#Separate out the controls
wetLabStarchCtrl<- wetLabStarch[which(wetLabStarch$Sample_ID == "control"),c(5,9:11)]
summary(wetLabStarchCtrl)
SD(wetLabStarchCtrl)

#Establish File to Hold stats
wetLabFile <- paste0("Data/OutputtedData/WSMDP_Wetlab_stats.txt")

#Write Starch Stats to File
cat("Starch Wetlab Data", file=wetLabFile, sep="\n", append=FALSE)
out <- capture.output(summary(wetLabStarchNoCtrl))
cat(out, file=wetLabFile, sep="\n", append=TRUE)
out <- capture.output(SD(wetLabStarchNoCtrl))
cat(out, file=wetLabFile, sep="\n", append=TRUE)
cat("Starch Wetlab Controls", file=wetLabFile, sep="\n", append=TRUE)
out <- capture.output(summary(wetLabStarchCtrl))
cat(out, file=wetLabFile, sep="\n", append=TRUE)
out <- capture.output(SD(wetLabStarchCtrl))
cat(out, file=wetLabFile, sep="\n", append=TRUE)

#########################
###Read in Wetlab Sugar Data ###
#########################

wetLabSugar <- read.csv("Data/RawData/WSMDP_Wetlab_Sugars.csv")

#set columns from char to numeric
cols <- c("Ave.Glucose","Ave.Fructose","Ave.Sucro","Ave.Total")
wetLabSugar[cols] <- lapply(wetLabSugar[cols], as.numeric)

#Separate out just useful data
wetLabSugarNoCtrl<- wetLabSugar[which(wetLabSugar$Variety_Num != "Control"),c(11,15:23,25:33)]
wetLabSugarNoCtrl<- wetLabSugarNoCtrl[which(wetLabSugarNoCtrl$Ave.Total != ""),]
Glucose <- c(wetLabSugarNoCtrl$Rep.A.Glucose, wetLabSugarNoCtrl$Rep.B.Glucose, wetLabSugarNoCtrl$Rep.C.Glucose )
Fructose <- c(wetLabSugarNoCtrl$Rep.A.Fructose, wetLabSugarNoCtrl$Rep.B.Fructose, wetLabSugarNoCtrl$Rep.C.Fructose )
Sucro <- c(wetLabSugarNoCtrl$Rep.A.Sucro, wetLabSugarNoCtrl$Rep.B.Sucro, wetLabSugarNoCtrl$Rep.C.Sucro )
Total <- c(wetLabSugarNoCtrl$Rep.A.Total, wetLabSugarNoCtrl$Rep.B.Total, wetLabSugarNoCtrl$Rep.C.Total )
wetLabSugarConj <- data.frame(Glucose,Fructose,Sucro,Total)

#Separate out just controls
wetLabSugarCtrl<- wetLabSugar[which(wetLabSugar$Variety_Num == "Control"),c(11,15:23,25:33)]
Glucose <- c(wetLabSugarCtrl$Rep.A.Glucose, wetLabSugarCtrl$Rep.B.Glucose, wetLabSugarCtrl$Rep.C.Glucose )
Fructose <- c(wetLabSugarCtrl$Rep.A.Fructose, wetLabSugarCtrl$Rep.B.Fructose, wetLabSugarCtrl$Rep.C.Fructose )
Sucro <- c(wetLabSugarCtrl$Rep.A.Sucro, wetLabSugarCtrl$Rep.B.Sucro, wetLabSugarCtrl$Rep.C.Sucro )
Total <- c(wetLabSugarCtrl$Rep.A.Total, wetLabSugarCtrl$Rep.B.Total, wetLabSugarCtrl$Rep.C.Total )
wetLabSugarControlConj <- data.frame(Glucose,Fructose,Sucro,Total)

summary(wetLabSugarControlConj)
SD(wetLabSugarControlConj)


#Write Sugar Stats to File
cat("Sugar Wetlab Data", file=wetLabFile, sep="\n", append=TRUE)
out <- capture.output(summary(wetLabSugarConj))
cat(out, file=wetLabFile, sep="\n", append=TRUE)
out <- capture.output(SD(wetLabSugarConj))
cat(out, file=wetLabFile, sep="\n", append=TRUE)
cat("Sugar Wetlab Controls", file=wetLabFile, sep="\n", append=TRUE)
out <- capture.output(summary(wetLabSugarControlConj))
cat(out, file=wetLabFile, sep="\n", append=TRUE)
out <- capture.output(SD(wetLabSugarControlConj))
cat(out, file=wetLabFile, sep="\n", append=TRUE)


WetLabSugarExpanded <- InfoCombination(wetLabSugar,SampleInfo,BookInfo)

#########################
###Read in the output from NIR ###
#########################

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
  # NIRDF$NIRBase<-gsub(NIRDF$NIR_ID, pattern = "-.*", replacement = "")
  OutDF$endo <- SampleInfoDF$endo[match(NIRDF$NIRBase, SampleInfoDF$NIRBase, nomatch = NA)]
  #eliminate SE markers
  OutDF$endo[which(OutDF$endo == "SE")] <- "field"
  OutDF$Inbred <- SampleInfoDF$Variety[match(NIRDF$NIRBase, SampleInfoDF$NIRBase, nomatch = NA)]
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

  #set everything to factors thats not carb data
  c1 <- which(colnames(OutDF)=="Starch")
  c2 <- which(colnames(OutDF)=="Total.Sugar")
  OutDF[-c(c1:c2)] <- lapply(OutDF[-c(c1:c2)], as.factor)
  return(OutDF)
}

HWSPsDF <- InfoCombination(HWSPsDF,SampleInfo,BookInfo)
LWSPWFsDF <- InfoCombination(LWSPWFsDF,SampleInfo,BookInfo)
LWSPNFsDF <- InfoCombination(LWSPNFsDF,SampleInfo,BookInfo)

#Eliminate irrelivant endosperm mutants from each df. We don't want to include sh2 samples that were predicted by the HWSP calibrated equations

HWSPs <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "su1se1"),]
LWSPWFs <- LWSPWFsDF[which(LWSPNFsDF$endo != "su1se1" & LWSPNFsDF$endo != "su1" & LWSPNFsDF$endo != ""),]
LWSPNFs <- LWSPNFsDF[which(LWSPNFsDF$endo == "sh2" | LWSPNFsDF$endo == "su1sh2-i"),]



# #lets look at the validation data that I predicted from kahtleen and Jareds data and output it for further analysis
# jared <- HWSPsDF[which(HWSPsDF$Year == "sc"),]
# write.csv(jared, file = "Data/OutputtedData/JaredsPredictedLineswHWSPeqn.csv")
# kathleen <- HWSPsDF[which(HWSPsDF$Year != "13"&HWSPsDF$Year != "14"&HWSPsDF$Year != "15"&HWSPsDF$Year != "sc"),]
# write.csv(kathleen,file = "Data/OutputtedData/KathleenPredictedLineswHWSPeqn.csv")
# 
########mash it all together ######
#All the predictions with the field lines
CarbInfoExpandedWFDF <- subset(rbind(HWSPs,LWSPWFs), !is.na(IsExperimental))
#all the predictions without the field lines
CarbInfoExpandedNFDF <- subset(rbind(HWSPs,LWSPNFs), !is.na(IsExperimental))

write.csv(file = "Data/OutputtedData/UncleanedInfoWFOutput.csv",CarbInfoExpandedWFDF)
write.csv(file = "Data/OutputtedData/UncleanedInfoNFOutput.csv",CarbInfoExpandedNFDF)

#visualize these data sets pre cleaning
CarbDataFrameVis(CarbInfoExpandedWFDF,"WithField_WithOutliers")
CarbDataFrameVis(CarbInfoExpandedNFDF,"NoField_WithOutliers")

#clean up the predictive data frames. Reassign or delete outliers
CleanedInfoWF <- CarbOutlierCleanup(CarbInfoExpandedWFDF,"WF",alpha = 0.05)
CleanedInfoNF <- CarbOutlierCleanup(CarbInfoExpandedNFDF,"NF",alpha = 0.05)

write.csv(file = "Data/OutputtedData/CleanedInfoWFOutput.csv",CleanedInfoWF)
write.csv(file = "Data/OutputtedData/CleanedInfoNFOutput.csv",CleanedInfoNF)

#revisualize the dataframes
CarbDataFrameVis(CleanedInfoWF,"WithField_Cleaned")
CarbDataFrameVis(CleanedInfoNF,"NoField_Cleaned")

#write the names of the varieties used to a csv file so we can find the corresponding GBS data
write.csv(file = "Data/OutputtedData/InbredsWithinWSMDPCarbDataNF.csv",unique(CleanedInfoNF[c("Inbred", "endo")]))
write.csv(file = "Data/OutputtedData/InbredsWithinWSMDPCarbDataWF.csv",unique(CleanedInfoWF[c("Inbred", "endo")]))


#clean up the predictive seperated data frames. Reassign or delete outliers
CleanedInfoHWSP <- CarbOutlierCleanup(HWSPs,"HWSP",alpha = 0.05)
CleanedInfoLWSPNF <- CarbOutlierCleanup(LWSPNFs,"LWSPNFs",alpha = 0.05)
CleanedInfoLWSPWF <- CarbOutlierCleanup(LWSPWFs,"LWSPWFs",alpha = 0.05)





#Take out the medians of each sample based on inbred and environment and visualize that too
CarbDFMedian <- CarbInfoExpandedWFDF %>%  
  group_by(Inbred, Envi) %>% 
  summarise(Starch = median(Starch)) 

CarbDFMedian$Total.Polysaccharides <- pull((CarbInfoExpandedWFDF %>%  
                                         group_by(Inbred, Envi) %>% 
                                         summarise(Total.Polysaccharides = median(Total.Polysaccharides)))[3] )

CarbDFMedian$WSP <- pull((CarbInfoExpandedWFDF %>%  
                       group_by(Inbred, Envi) %>% 
                       summarise(WSP = median(WSP)))[3] )
CarbDFMedian$Glucose <- pull((CarbInfoExpandedWFDF %>%  
                           group_by(Inbred, Envi) %>% 
                           summarise(Glucose = median(Glucose)))[3] )
CarbDFMedian$Fructose <- pull((CarbInfoExpandedWFDF %>%  
                            group_by(Inbred, Envi) %>% 
                            summarise(Fructose = median(Fructose)))[3] )
CarbDFMedian$Sucrose <- pull((CarbInfoExpandedWFDF %>%  
                           group_by(Inbred, Envi) %>% 
                           summarise(Sucrose = median(Sucrose)))[3] )
CarbDFMedian$Total.Sugar <- pull((CarbInfoExpandedWFDF %>%  
                               group_by(Inbred, Envi) %>% 
                               summarise(Total.Sugar = median(Total.Sugar)))[3] )
CarbDFMedian$endo <- pull((CarbInfoExpandedWFDF %>%  
                                    group_by(Inbred, Envi) %>% 
                                    summarise(Total.Sugar = last(endo)))[3] )

# CarbDFMedian[c(3:9)] <- lapply(CarbDFMedian[c(3:9)], as.numeric)

CarbDataFrameVis(CarbDFMedian,"WithField_WithOutliers_Inbred-enviMedian")



##Now look at the median for just the inbreds
CarbDFMedian <- CarbInfoExpandedWFDF %>%  
  group_by(Inbred) %>% 
  summarise(Starch = median(Starch)) 

CarbDFMedian$Total.Polysaccharides <- pull((CarbInfoExpandedWFDF %>%  
                                              group_by(Inbred) %>% 
                                              summarise(Total.Polysaccharides = median(Total.Polysaccharides)))[2] )

CarbDFMedian$WSP <- pull((CarbInfoExpandedWFDF %>%  
                            group_by(Inbred) %>% 
                            summarise(WSP = median(WSP)))[2] )
CarbDFMedian$Glucose <- pull((CarbInfoExpandedWFDF %>%  
                                group_by(Inbred) %>% 
                                summarise(Glucose = median(Glucose)))[2] )
CarbDFMedian$Fructose <- pull((CarbInfoExpandedWFDF %>%  
                                 group_by(Inbred) %>% 
                                 summarise(Fructose = median(Fructose)))[2] )
CarbDFMedian$Sucrose <- pull((CarbInfoExpandedWFDF %>%  
                                group_by(Inbred) %>% 
                                summarise(Sucrose = median(Sucrose)))[2] )
CarbDFMedian$Total.Sugar <- pull((CarbInfoExpandedWFDF %>%  
                                    group_by(Inbred) %>% 
                                    summarise(Total.Sugar = median(Total.Sugar)))[2] )
CarbDFMedian$endo <- pull((CarbInfoExpandedWFDF %>%  
                             group_by(Inbred) %>% 
                             summarise(Total.Sugar = last(endo)))[2] )

# CarbDFMedian[c(3:9)] <- lapply(CarbDFMedian[c(3:9)], as.numeric)

#need to do this manually since no envi data
# CarbDataFrameVis(CarbDFMedian,"WithField_WithOutliers_Inbred-enviMedian")


#######################
### Summarize the samples within CleanedInfoNF
summarizeDF <- function(DF){
print(summary(DF))
print(length(unique(DF$Inbred)))
OnlyUnique <- DF %>% distinct(Inbred, .keep_all = TRUE)
print(summary(OnlyUnique))
}
summarizeDF(CleanedInfoWF)
summarizeDF(CleanedInfoNF)
checkfile <- "Data/OutputtedData/CheckDistributions.csv"
checkfile2 <- "Data/OutputtedData/CheckDistributionsByEnvi.csv"

CleanedInfoWFCheck1 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check1"),c(5:11,15)]
carbs <- matrix(c("",colnames(CleanedInfoWF)[5:11]),nrow = 1, ncol = 8)
write.csv(carbs,checkfile,append = TRUE, row.names = FALSE)
write.csv("Check 1, We05407",checkfile,append = TRUE, row.names = FALSE)
write.csv(summary(CleanedInfoWFCheck1), checkfile,append = TRUE)
write.csv(carbs,checkfile,append = TRUE, row.names = FALSE)
write.csv("Standard Deviation",checkfile,append = TRUE, row.names = FALSE)
write.csv(SD(CleanedInfoWFCheck1), checkfile,append = TRUE)
# hist(CleanedInfoWFCheck2$Starch)

Check1 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check1"),]
out1 <- Check1 %>%
  group_by(Envi) %>%
  summarise_at(vars(5:11), list(name = mean), na.rm = TRUE)
carbs <- matrix(c("Envi",colnames(CleanedInfoWF)[5:11]),nrow = 1, ncol = 8)
write.csv(carbs,checkfile2,append = TRUE, row.names = FALSE)
write.csv("Check 1, We05407",checkfile2,append = TRUE, row.names = FALSE)
write.csv(out1, checkfile2,append = TRUE, row.names = FALSE)


CleanedInfoWFCheck2 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check2"),c(5:11,15)]
write.csv("Check 2, W5579",checkfile,append = TRUE, row.names = FALSE)
write.csv(summary(CleanedInfoWFCheck2), checkfile,append = TRUE)
write.csv("Standard Deviation",checkfile,append = TRUE, row.names = FALSE)
          write.csv(SD(CleanedInfoWFCheck2), checkfile,append = TRUE)
Check2 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check2"),]
out2 <- Check2 %>%
  group_by(Envi) %>%
  summarise_at(vars(5:11), list(name = mean), na.rm = TRUE)
write.csv("Check 2, W5579",checkfile2,append = TRUE, row.names = FALSE)
write.csv(out2, checkfile2,append = TRUE, row.names = FALSE)
          
          
CleanedInfoWFCheck3 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check3"),c(5:11,15)]
write.csv("Check 3, Ia5125",checkfile,append = TRUE, row.names = FALSE)
write.csv(summary(CleanedInfoWFCheck3), checkfile,append = TRUE)
write.csv("Standard Deviation",checkfile,append = TRUE, row.names = FALSE)
          write.csv(SD(CleanedInfoWFCheck3), checkfile,append = TRUE)
Check3 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check3"),]
out3 <- Check3 %>%
  group_by(Envi) %>%
  summarise_at(vars(5:11), list(name = mean))
write.csv("Check 3, Ia5125",checkfile2,append = TRUE, row.names = FALSE)
write.csv(out3, checkfile2,append = TRUE, row.names = FALSE)

CleanedInfoWFCheck4 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check4"),c(5:11,15)]
write.csv("Check 4, IL125b",checkfile,append = TRUE, row.names = FALSE)
write.csv(summary(CleanedInfoWFCheck4), checkfile,append = TRUE)
write.csv("Standard Deviation",checkfile,append = TRUE, row.names = FALSE)
          write.csv(SD(CleanedInfoWFCheck4), checkfile,append = TRUE)
Check4 <- CleanedInfoWF[which(CleanedInfoWF$BookInbred == "Check4"),]
out4 <- Check4 %>%
  group_by(Envi) %>%
  summarise_at(vars(5:11), list(name = mean))
write.csv("Check 4, IL125b",checkfile2,append = TRUE, row.names = FALSE)
write.csv(out4, checkfile2,append = TRUE, row.names = FALSE)
###Validate that the NIR Equation is good###
#########################

#######Equation Validation!##############
#Now I have a variable that has all the projected values, for the values used to calibrate the equations. WHat are the statistics on that?
#set carb to a number 1:7. carb <- c(Fructose, Glucose, Sucrose, Total Sugar, Starch, Total Polysaccharide, WSP)
#this functions needs to have the data frame set up so the carbs being compared are directly next to eachother
EqnStats <- function(DF){

  #establish dataframe used to record stats
  Out <- data.frame(Carb = c("Starch", "Total Polysaccharide", "WSP", "Glucose","Fructose", "Sucrose", "Total Sugar"),
                    RMSEP = rep(NA,7),
                    bias = rep(NA,7),
                    SEE = rep(NA,7),
                    slope = rep(NA,7),
                    intercept = rep(NA,7),
                    R2 = rep(NA,7))
  
  #look at each trait in turn

  dfpos <- c(2,4,6,8,10,12,14)
  for(carb in 1:7){
    #Calculated the RMSEP
    Out$RMSEP[carb] <- sqrt(sum((DF[,dfpos[carb]]- DF[,dfpos[carb]-1])^2, na.rm = TRUE)/dim(DF)[1])
    #Calculate the Bias
    Out$bias[carb] <- mean(DF[,dfpos[carb]], na.rm = TRUE) - mean(DF[,dfpos[carb]-1],na.rm = TRUE)
    #Calculated the SEE
    Out$SEE[carb] <- sqrt((dim(DF)[1]/(dim(DF)[1]-1))*(Out$RMSEP[carb]^2-Out$bias[carb]^2))
    # Out$Carb[carb] <- colnames(DF[dfpos[carb]-1])
  }
  return(Out)
}


#Visualize these
R2Vis <- function(DF, label, Out){
  dfpos <- c(2,4,6,8,10,12,14)
  alpha <- LETTERS[1:8]
  Carb = c("Starch", "Total Polysaccharide", "WSP", "Glucose","Fructose", "Sucrose", "Total Sugar")
  carbFileName <- paste("Figures/wsmdp2021_",label,"_AllCarb_NIREqn_Prediction_vis.pdf", sep = "")
  pdf(carbFileName, height= 14)
  par(mfrow=c(4,2))
  for(i in 1:7){
    carbCompare  <- lm(DF[,dfpos[i]]~ DF[,dfpos[i]-1])
    print(summary(carbCompare))
    rsqua <- round(summary(carbCompare)$r.squared,2)
    plot(DF[,dfpos[i]]~ DF[,dfpos[i]-1],
         pch = 16,
         xlab = "Assay Results (%)",
         ylab = "NIR Prediction (%)",
         # main = expression(Carb[i]*'r'^2*" ="*!!rsqua))
         main = paste0(alpha[i],": ",Carb[i]))
    abline(coefficients(carbCompare), lwd = 2, lty = 2, col = "red")
   
    legend("bottomright",paste0("r^2 =",rsqua),  
           box.lty=0)
    # text(15,max(Prediction[,i])-5,labels = paste("r^2 =",trunc(rsqua*10^3)/10^3))

    Out$slope[i] <- trunc(10^3*summary(carbCompare)$coefficients[2])/10^3
    Out$intercept[i] <- trunc(10^3*summary(carbCompare)$coefficients[1])/10^3
    Out$R2[i] <- rsqua
    
  }
  dev.off()
  return(Out)
}


######Visualize the validation. Look at all the samples that were wetlabbed and used to calibrate the equation
wetlab <- read.csv("Data/WSMDP_Wetlab_StarchSugarData_FormatedForWinISI_WithR.csv")
wetlabDF <- wetlab[,c(1,6:8,3,2,4,5)]

WLLabels <- c("Samples","Starch_WL", "Total Polysaccharide_WL", "WSP_WL", "Glucose_WL","Fructose_WL", "Sucrose_WL", "Total Sugar_WL")
colnames(wetlabDF) <- WLLabels

#merge the actual and the predicted traits
WFWL <- merge(CleanedInfoWF, wetlabDF, by = "Samples")
NFWL <- merge(CleanedInfoNF, wetlabDF, by = "Samples")

#this pulls out just the traits and eliminates extraneous info
alternatingorder <- c(1,5,28,6,29,7,30,8,31,9,32,10,33,11,34)

#this makes a data frame where the traits alternate. Predicted starch, actual (wetlab) starch, Predicted totalpoly, actual totalpoly, etc etc
WFWLdf <- WFWL[,alternatingorder]
NFWLdf <- NFWL[,alternatingorder]

#calculate some of the statistics about this dataframe
EqnStatisticOut <- EqnStats(WFWLdf[,2:15])
#read the output from the eqnstatistics and plug it into the R2Vis function. string is name for saving the files
WFWLdfEqnStatsR <- R2Vis(WFWLdf[,2:15], "CleanedWSPeqnWF_PredVsWetlab_for_Calibration_Samples", EqnStatisticOut)
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


# ####With Jared WetlabDATA
# HWSPsJ <- HWSPsDF[which(HWSPsDF$endo == "su1" | HWSPsDF$endo == "se" | is.na(HWSPsDF$endo)),]
# alternatingorder2 <- c(1,5,27,6,28,7,29,8,30,9,31,10,32,11,33)
# CarbInfoExpandedWFJDF <- rbind(HWSPsJ,LWSPWFs)
# CarbInfoExpandedNFJDF <- rbind(HWSPsJ,LWSPNFs)
# WFValWLJ <- merge(CarbInfoExpandedWFJDF, valwetlabDF, by = "Samples")
# NFValWLJ <- merge(CarbInfoExpandedNFJDF, valwetlabDF, by = "Samples")
# WFValWLJdf <- WFValWLJ[19:30,alternatingorder2]
# NFValWLJdf <- NFValWLJ[19:30,alternatingorder2]
# WFValWLdfJEqnStatsR <- R2Vis(WFValWLJdf[,2:15], "UnCleanedWSPeqnWF_PredVsWetlab_ValidationSubset_WJared", EqnStats(WFValWLdf[,2:15]))
# NFValWLdfJEqnStatsR <- R2Vis(NFValWLJdf[,2:15], "UnCleanedWSPeqnNF_PredVsWetlab_ValidationSubset_WJared", EqnStats(NFValWLdf[,2:15]))
# 
# write.csv(WFValWLdfJEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationWF_WJared.csv")
# write.csv(NFValWLdfJEqnStatsR, "Data/OutputtedData/EqnFitStatisticsValidationNF_WJared.csv")

#########################
###Linear Model Analysis###
#########################

linearmodel <- function(SampleDFtoModel,TitleAddendum, endoCheck = FALSE){
  ptm <- proc.time()
  alpha <- LETTERS[1:8]
  #For Debugging
  # SampleDFtoModel <- CleanedInfoNF
  # TitleAddendum <- "CleanedOutliersNF"
  #Establish a file to output the statistical analysis results to
  statsfile <- paste0("Data/OutputtedData/WSMDP_CarbPheno_stats_",TitleAddendum,".txt")
  #write a header in that stats file
  cat(paste0("Phenotypic Statistics for the Carbohydrates of ",TitleAddendum," Through Mixed Linear Modeling."), file=statsfile, sep="\n", append=FALSE)

  #make a dataframe to hold the BLUPs after the models are done
  blupHolder <- data.frame("Inbred" = sort(unique(SampleDFtoModel$BookInbred)))
  #make a dataframe to hold the variances due to different factors
  VarDF <- data.frame("Carb" = colnames(SampleDFtoModel)[5:11],"Inbred" = rep(NA,7),"Envi" = rep(NA,7),"Inbred:Envi" = rep(NA,7),  "superblock" = rep(NA,7), "Col" = rep(NA,7),"Row" = rep(NA,7),"block" = rep(NA,7),"Residuals"= rep(NA,7))
  
  
  #figure out where the traits start and stop
  c1 <- which(colnames(SampleDFtoModel)=="Starch")
  c2 <- which(colnames(SampleDFtoModel)=="Total.Sugar")
  #write those traits to a vector
  carbs <- colnames(SampleDFtoModel)[c1:c2]
  
  ######Iteratite though all the traits#####
  pdf(paste0("Figures/WSMDP_LinearModel_assumptions_",TitleAddendum,"_AllCarbs.pdf"), width = 6, height = 12)
  par(mfrow=c(7,3))
  for(j in 1:length(carbs)){
    print(paste0("Starting ",carbs[j]," Model Generation"))
    #output the statistics about this trait to the stats file
    summary <- summary(SampleDFtoModel[,c1+j-1], )
    out <- capture.output(summary)
    cat(out, file=statsfile, sep="\n", append=TRUE)

    #####Establish Model
    modelpastecheck<-  paste0(carbs[j], " ~ (1|BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")

    # Use check if there enough checks in the model set. Unique numbers if no check are NA and 99
    if(length(unique(SampleDFtoModel$Check))>2){
    modelpastecheck<-  paste0(carbs[j], " ~ Check + (1|BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")
    }
    
  #Output the model to the stats file
  cat(modelpastecheck, file=statsfile, sep="\n", append=TRUE)

  #Run the Model
  #get rid of the one sample that is throwing an error because its NA
  SampleDFtoModel_noTraitNAs<-subset(SampleDFtoModel, !is.na(superblock))   
  # SampleDFtoModel_noTraitNAs <- SampleDFtoModel_noTraitNAs[-which(is.na(SampleDFtoModel_noTraitNAs[carbs[j]])),]
  #plug in the model
  print("Run the Big Model")
  modelExpansive <- lme4::lmer(modelpastecheck,
                data=SampleDFtoModel, REML = TRUE)
  #Step through the model to find the shortened, better fitting version
  print(paste0("Step Through ",carbs[j]," Model Generation"))
  x <- stepcAIC(modelExpansive,direction="backward")
  #Extract that shorter model
  model <- x$finalModel
  
  print(paste0("Read Out ",carbs[j]," BLUPS from selected model"))
  #get the random effects of that model. This is deviation from the average
  RandomEffects <- ranef(model)
  #get mean of the trait
  average <- mean(pull(SampleDFtoModel[,c1+j-1]),na.rm = T)
  #add the average to get the predicted trait value, not deviation of the book inbreds
  BLUPofTrait <- RandomEffects$BookInbred + average 
  #hold the blups
  # tempBlup <- data.frame("Inbred" = rownames(RandomEffects$BookInbred), "BLUP" = RandomEffects$BookInbred)
  tempBlup <- data.frame("Inbred" = rownames(RandomEffects$BookInbred), "BLUP" = BLUPofTrait)
  blupHolder <- merge(blupHolder, tempBlup, by = "Inbred", all = TRUE)
  colnames(blupHolder)[j+1] <- paste0(carbs[j],".BLUP")
  

  
  ###### Check Assumptions#######
  print(paste0("Check  ",carbs[j]," Model Assumption"))
  # Model Fit with REML
  plot(fitted(model), residuals(model), pch=19, col="dark blue", ylab="Residuals", xlab="Predicted")
  abline(h=0,col="red", lwd=1, lty=1)
  # histogram of residuals
  hist(residuals(model),main=paste0(alpha[j], ": ", carbs[j]," residuals"),freq=F, xlab="Residuals", ylab= "Freq", col="palegreen", col.main="darkblue")
  x=seq(-5e-15,9e-15,5e-15)
  curve(dnorm(x,mean(residuals(model)),sd(residuals(model))),add=T,lwd=2, col="red", lty=1)
  # qq plot
  qqnorm(residuals(model), pch=19, col="dark blue", col.lines="red", xlab="Pred quantiles", ylab="Obs quantiles", main = "")


  # Summary of random effects
  summary <- summary(model, correlation=FALSE)
  out <- capture.output(summary)
  cat(out, file=statsfile, sep="\n", append=TRUE)
  
  # Write out residuals from ANOVA
  write.table(resid(model), paste0("Data/OutputtedData/WSMDP_LinearModel_residuals_",TitleAddendum,"_", carbs[j], ".csv"), col.names=F, row.names=F, sep=",")

  #######calculate model R2
  print(paste0("Calculate ",carbs[j]," R2"))
  r2 <- MuMIn::r.squaredGLMM(model)
  out <- capture.output(r2[2])
  cat(paste("The r^2 value of this model for", carbs[j] ,"is", out), file= statsfile, sep="\n", append=TRUE)
    
  
  ######## Calculate heritability
  model_variances <- as.data.frame(VarCorr(model))
  print(paste0("Calculate ",carbs[j]," Heritability"))
  #Where are the geno, gxe, and residual stats located in the dataframe?
  geno <- which(model_variances$grp == "BookInbred")
  gxe <- which(model_variances$grp == "BookInbred:Envi")
  resid <- which(model_variances$grp == "Residual")

  #broad sence heritability calculated as variance of genotype/ (variance of geno + var of gxe / num of enviornments + var of error/num of replicates and environ)
  h2 <- model_variances$vcov[geno]/(model_variances$vcov[geno]+(model_variances$vcov[gxe]/4)+(model_variances$vcov[resid]/8))
  out <- capture.output(h2)
  cat(paste("The heritability of", carbs[j] ,"is", out), file= statsfile, sep="\n", append=TRUE)
  
  #### Shapiro-Wilk test normality test 
  ####p-val < alpha : reject H~0~ that the residuals are normally distributed.
  out <- shapiro.test(residuals(model))
  print(out)
  cat(paste("The Shapiro-Wilk test normality test results of", carbs[j] ,"is", out, ". p-val < alpha : reject H~0~ that the residuals are normally distributed."), file= statsfile, sep="\n", append=TRUE)
  
  # ###### Levene's test 
  # ###Will provide a p-value to unequal variance testing.p-val < alpha : reject H~0~ that the variances are equal.
  # out <- leveneTest(carbs[j]~endo,data=SampleDFtoModel)
  # print(out)
  # cat(paste("The Levene's test results of", carbs[j] ," and endosperm mutant is", out, ". Will provide a p-value to unequal variance testing.p-val < alpha : reject H~0~ that the variances are equal."), file= statsfile, sep="\n", append=TRUE)
  # 
  #Only Do the parts with the endosperm component if endoCheck is set to true since it takes a long time
  if(endoCheck == TRUE){
  ######## Compare model with one with Endosperm #####
  cat("Model Comparison, with and without Endosperm Term", file=statsfile, sep="\n", append=TRUE)
  #####Modeling
  print(paste0("Check ",carbs[j]," Model with Endosperm components"))
  modelpasteendo<-  paste0(carbs[j], " ~ (1|endo/BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")
  # Use check if there enough checks in the model set
  if(length(unique(SampleDFtoModel$Check))>1){
  modelpasteendo<-  paste0(carbs[j], " ~ Check + (1|endo/BookInbred) + (1|BookInbred:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col)")
  }
  cat(modelpasteendo, file=statsfile, sep="\n", append=TRUE)

  #plug in the model
  modelendoExpansive <- lmer(modelpasteendo,
                         data=SampleDFtoModel_noTraitNAs, REML = TRUE)
  #Step through the model to find the shortened, better fitting version
  print(paste0("Step through ",carbs[j]," Model with Endosperm components"))
  x <- stepcAIC(modelExpansive,direction="backward")
  #Extract that shorter model
  modelendo <- x$finalModel
  
  # Summary of random effects
  summary <- summary(modelendo, correlation=FALSE)
  out <- capture.output(summary)
  cat(out, file=statsfile, sep="\n", append=TRUE)
  }
  
  
  
  print(paste0("Find what variances arribute to what in ",carbs[j]))
  formula1 <- paste0(colnames(SampleDFtoModel)[j+4],"~ BookInbred*Envi  + superblock%in%Envi + superblock%in%block%in%Envi+ Row%in%Envi+ Col%in%Envi")
  fit1 <- lm(formula1,data=SampleDFtoModel)
  AIC1 <- extractAIC(fit1)
  AIC1
  summary(fit1)$r.square
  print(anova(fit1))
  out <- anova(fit1)
  SStotal <- sum(out$`Sum Sq`)
  for(i in 1:8){
    #variance explained is caluclated by the sum of squares divided by the sum of squares total
    VarDF[j,i+1] <- (out$`Sum Sq`[i]/SStotal)}
  
  }
  dev.off()

  #######Graph the different variances explained by different factors######
  VarDFMelt <- reshape2::melt(VarDF)
  pdf(paste("Figures/WSMDP_AllNIRPred_MixedEqn_PercentVarianceExplainedby_Factors_",TitleAddendum,".pdf",sep=""), width = 10, height = 5)
  print(barchart(~value|variable, group = factor(Carb), data= VarDFMelt,reverse.rows = FALSE,main = "Percent Phenotypic Variance Explained",layout = c(8,1),
           key = simpleKey(text = colnames(SampleDFtoModel)[5:11],
                           rectangles = TRUE, points = FALSE, space = "right")))
  dev.off()
  
  #######Graph the different variances explained by different factors with GGPlot ######
  pdf(paste("Figures/WSMDP_AllNIRPred_MixedEqn_PercentVarianceExplainedby_Factors_",TitleAddendum,"_ggplot.pdf",sep=""), width = 10, height = 5)
  p <- ggplot(VarDFMelt, aes(y = Carb, value)) +
    geom_bar(aes(fill = variable),stat = "identity",  position = position_stack(reverse = TRUE)) + 
    theme(legend.position = "top",text = element_text(size=20))
  print(p)
  dev.off()
  
  #Write the variances to a file
  write.table(VarDF, file = paste0("Data/OutputtedData/WSMDP_CarbPheno_Anova_Variances_",TitleAddendum,".txt"), col.names=T, row.names=F, sep=",")
  
  
  
  # Write out BLUPs for Genotypes
  write.table(blupHolder, file=paste0("Data/OutputtedData/WSMDP_CarbPheno_InbredBLUPS_",TitleAddendum,".txt"), col.names=T, row.names=F, sep=",")
  beep(3)
  print(proc.time() - ptm)
  return(blupHolder)
}

# WFBlups <- linearmodel(CleanedInfoWF,"CleanedOutliersWF")
# NFBlups <- linearmodel(CleanedInfoNF,"CleanedOutliersNF",endoCheck = FALSE)
JustSu1 <- HWSPsDF[which(HWSPsDF$endo == "su1" ),]

JustSu1Blups <- linearmodel(JustSu1, "Uncleanedsu1only")
Justsh2Blups <- linearmodel(LWSPNFs, "Uncleanedsh2nionly")

# 
# HWSPBlups <- linearmodel(CleanedInfoHWSP,"CleanedOutliersHWSP")
# These dont work. Fix. Not enough of a certain factor type?
# LWSPNFBlups <- linearmodel(CleanedInfoLWSPNF,"CleanedOutliersLWSPNF")
# LWSPWFBlups <- linearmodel(CleanedInfoLWSPWF,"CleanedOutliersLWSPWF")





#########################
###SNP Relate###
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
# LDMatFile <- paste0("Figures/",Sys.Date(),"SNPRelate_LDMatrix.jpeg")
# jpeg(LDMatFile, width = 1000, height = 1000)
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

#########only include the samples that are in genoinfo########
#this is a round about way of looking at only the genotype info I have for the things planted in this trial
#eliminate inbreds tested in other experiments
genoinfo2014only <- genoinfo[which(genoinfo$Planting20142015 == 1),]

#look at what inbreds have the same genoname in both the seq file and the data summary
matching <- data.frame("GenoName" = unique(PCA$sample.id[na.omit(match(genoinfo2014only$GenoName, holdtrunc))]))

#hold the names of the inbreds I have seq info for
holdnocol <- unique(holdtrunc[na.omit(match(genoinfo2014only$GenoName, holdtrunc))])

#combine the geno name with the endosperm mutant type from the genoinfo file
matching$endo <- genoinfo2014only$endo[match(holdnocol,genoinfo2014only$GenoName)]

#combine the geno name with the endosperm mutant type from the genoinfo file
matching$Program <- genoinfo2014only$Program[match(holdnocol,genoinfo2014only$GenoName)]
#combine the geno name with the endosperm mutant type from the genoinfo file
matching$Region <- genoinfo2014only$Region[match(holdnocol,genoinfo2014only$GenoName)]


#conduct pca again with only those samples
PCA <- snpgdsPCA(PreLDGeno, sample.id = matching$GenoName, snp.id = PostLDGeno.id)
pc.perc <- PCA$varprop*100
head(round(pc.perc,2))


ptm <- proc.time()
#visualize the PCA with ld pruned snps and only the tested samples
PCAFigureCreation(PCA,pc.perc,matching,infilename,"endo")
PCAFigureCreation(PCA,pc.perc,matching,infilename,"Program")
PCAFigureCreation(PCA,pc.perc,matching,infilename,"Region")
proc.time() - ptm
#########################
### Close Snp Relate ###
#########################
outfilename <- "WSMDP_SeqG"
outfile <- paste("Data/RawData/",outfilename ,sep="")
snpgdsGDS2PED(PreLDGeno, outfile, snp.id = PostLDGeno.id)
snpgdsGDS2PED(PreLDGeno, paste0("Data/RawData/",outfilename,"_NoLD_InbredsPruned"), sample.id = matching$GenoName)
snpgdsClose(PreLDGeno)



#########################
### GWAS Time ###
#########################

#read in genetic info post MAF
Seq <- "SeqG"
hmppath <- paste0("Data/RawData/WSMDP_",Seq,".hmp.txt")


# hmppath <- "Data/RawData/WSMDP_SCMV_SeqB.hmp.txt"
# hmppath <- "Data/RawData/WSMDP_SeqC.hmp.txt"
geno <- fread(hmppath,skip = "rs#")
str(geno)
colnames(geno)<-gsub(colnames(geno), pattern = ":.*", replacement = "")
str(geno)
# geno <- read.csv("Data/RawData/SeqG_numericFormat.csv")

######## Read in blups from file########
#this way you don't have to rerun the blup generation linear model function every time. that takes a long time

BlupGenoCleanup <- function(BlupDFName){
BlupDF <- read.csv(file=paste0("Data/OutputtedData/WSMDP_CarbPheno_InbredBLUPS_CleanedOutliers",BlupDFName,".txt"))#, col.names=T, row.names=F)
  
BlupDFGeno <- merge(BlupDF,genoinfo, by = "Inbred")
c1 <- which(colnames(BlupDFGeno)=="Starch.BLUP")
c2 <- which(colnames(BlupDFGeno)=="Total.Sugar.BLUP")
c3 <- which(colnames(BlupDFGeno)=="GenoName")
c4 <- which(colnames(BlupDFGeno)=="endo")

#remove rows with no geno information and remove rows that are duplicates. 
#TODO:Check which duplicates are being deleted and if its a good choice
BlupDFGeno <- BlupDFGeno[-which(BlupDFGeno$GenoName == ""),]
# BlupDFGeno <- BlupDFGeno[-which(BlupDFGeno$GenoName == "0"),]
# BlupDFGeno <- BlupDFGeno[-which(is.na(BlupDFGeno$GenoName)),]
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
# HWSPBlupsGenoJustPheno <- BlupGenoCleanup("HWSP")
# LWSPNFBlupsGenoJustPheno <- BlupGenoCleanup("LWSPNF")
# LWSPWFBlupsGenoJustPheno <- BlupGenoCleanup("LWSPWF")

#Establish what and where the traits are for the gwas
c1 <- which(colnames(NFBlupsGenoJustPheno)=="Total.Polysaccharides.BLUP")
c2 <- which(colnames(NFBlupsGenoJustPheno)=="WSP.BLUP")
c1 <- which(colnames(NFBlupsGenoJustPheno)=="Starch.BLUP")
# c1 <- which(colnames(WFBlupsGenoJustPheno)=="Glucose.BLUP")
c2 <- which(colnames(NFBlupsGenoJustPheno)=="Total.Sugar.BLUP")
blups <- colnames(NFBlupsGenoJustPheno[c1:c2])



gwasPlotsNo <- vector('list', 7)
gwasPlotsEndo <- vector('list', 7)
count = 1
# Start the clock!
ptm <- proc.time()
#itterate through all the traits.
for(blup in blups){
  
  # gwasPlotsNo[[count]] <-  GWASPolyRunner(WFBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_",Sys.Date()),Seq,"WFBLUP")
  # gwasPlotsEndo[[count]] <- GWASPolyRunner(WFBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_",Sys.Date()),Seq,"WFBLUP","endo","factor")
  
  gwasPlotsNo[[count]] <- GWASPolyRunner(NFBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"NFBLUP",NULL,NULL,Thresh = "M.eff")
  gwasPlotsEndo[[count]] <- GWASPolyRunner(NFBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"NFBLUP","endo","factor",Thresh = "M.eff")
  
  # GWASPolyRunner(HWSPBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"WFBLUP")
  # GWASPolyRunner(HWSPBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"WFBLUP","endo","factor")
  
  # GWASPolyRunner(LWSPNFBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPNFBlups")
  # GWASPolyRunner(LWSPNFBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPNFBlups","endo","factor")
  # 
  # GWASPolyRunner(LWSPWFBlupsGenoJustPheno[,1:8],geno,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPWFBlups")
  # GWASPolyRunner(LWSPWFBlupsGenoJustPheno,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"LWSPWFBlups","endo","factor")
  
  count <- count + 1
}

beep(2)
#Stop the clock!
proc.time() - ptm

removelegend <- function(p){
  g <- p +  theme(legend.position = "none")
  return(g)
}

gwasPlotsNope <- vector('list', 7)
gwasPlotsEndope <- vector('list', 7)
gwasPlotsNope <- lapply(gwasPlotsNo, removelegend )
gwasPlotsEndope <- lapply(gwasPlotsEndo, removelegend )

# pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoFixedEffect_FDRThresh_2021-10-20_NFBLUP_ALLTraitsGWASCumulative_GeneralModel.pdf",width = 10, height = 10)
# pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoFixedEffect_FDRThresh_2021-12-09_NFBLUP_ALLTraitsGWASCumulative_AdditiveModel.pdf",width = 10, height = 10)
pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoFixedEffect_FDRThresh_2021-12-16_NFBLUP_ALLTraitsGWASCumulative_AdditiveModel.pdf",width = 10, height = 10)
figure <- ggarrange(gwasPlotsNope[[1]], gwasPlotsNope[[2]], gwasPlotsNope[[3]], gwasPlotsNope[[4]], gwasPlotsNope[[5]], gwasPlotsNope[[6]], gwasPlotsNope[[7]],
                    labels = c("A", "B","C","D","E","F","G"),
                    ncol = 2, nrow = 4,
                    font.label = list(size = 18))

print(figure)
dev.off()

# pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_EndoFixedEffect_FDRThresh_2021-10-20_NFBLUP_ALLTraitsGWASCumulative_GeneralModel.pdf",width = 10, height = 10)
# pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_EndoFixedEffect_FDRThresh_2021-12-09_NFBLUP_ALLTraitsGWASCumulative_AdditiveModel.pdf",width = 10, height = 10)
pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_EndoFixedEffect_FDRThresh_2021-12-16_NFBLUP_ALLTraitsGWASCumulative_AdditiveModel.pdf",width = 10, height = 10)
figure <- ggarrange(gwasPlotsEndope[[1]], gwasPlotsEndope[[2]], gwasPlotsEndope[[3]], gwasPlotsEndope[[4]], gwasPlotsEndope[[5]], gwasPlotsEndope[[6]], gwasPlotsEndope[[7]],
                    labels = c("A", "B","C","D","E","F","G"),
                    ncol = 2, nrow = 4,
                    font.label = list(size = 18))

print(figure)
dev.off()

# pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoandEndoFixedEffect_FDRThresh_2021-10-20_NFBLUP_TotalSugar_GeneralModel.pdf",width = 10, height = 10)
# pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoandEndoFixedEffect_FDRThresh_2021-12-09_NFBLUP_TotalSugar_AdditiveModel.pdf",width = 10, height = 10)
pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoandEndoFixedEffect_FDRThresh_2021-12-16_NFBLUP_TotalSugar_AdditiveModel.pdf",width = 10, height = 10)
figure <- ggarrange(gwasPlotsNope[[7]],gwasPlotsEndope[[7]],
                    labels = c("A","B"),
                    ncol = 1, nrow = 2,
                    font.label = list(size = 18))

print(figure)
dev.off()

#########################
###plot the gwas results###
#########################
#read in QTL results from the files made by GWASPolyRunner and GWASPolyVis
# GWASPolyRunVersion <- paste0("EndoFixedEffect_FDRThresh_",Sys.Date())

FullGWASVisualize <- function(GWASPolyRunVersion, DataSet = "NFBLUP", Thresh = "FDR"){
  # for debugging
  # Thresh = "FDR"
  # DataSet = "NFBLUP"
  # # GWASPolyRunVersion <- "EndoFixedEffect_FDRThresh_2021-10-20"
  # GWASPolyRunVersion <- "EndoFixedEffect_FDRThresh_2021-12-09"
  # Seq <- "SeqG"
  # hmppath <- paste0("Data/RawData/WSMDP_",Seq,".hmp.txt")
  # blups <- colnames(NFBlupsGenoJustPheno[2:8])
  
  #establish DF to hold the file readins 
  QTLList <- list()


#Put these files into a dataframe
#This is the originally created filename
for(i in 1:7){
  file1 <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,DataSet,GWASPolyRunVersion,"_",blups[i],"_",Thresh,"_FitQTLLoci.csv", sep = "")
  file2 <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,"_",DataSet,"_",GWASPolyRunVersion,"_",blups[i],"_",Thresh,"_SignificantLoci.csv", sep = "")
  #make sure theres actually something in the file: As in, there are any QTLs and you don't try to open an empty file
  if (file.exists(file1) && file.size(file1)>15){ 
    merge1 <- read.csv(file1)
    merge2 <- read.csv(file2)
    if(dim(merge1)[1] != dim(merge2)[1]){
      stop("Something is wrong with the QTL files")
    }
    merge <- merge(merge1, merge2)
  QTLList[[i]] <- merge}}
#bind the QTLs into one dataframe
QTLDF <- bind_rows(QTLList)
#pull out only the general QTL
# QTLDFGen <- QTLDF[which(QTLDF$Model == "general"),]
QTLDFGen <- QTLDF[which(QTLDF$Model == "additive"),]
QTLDFGen$Trait <- paste0(QTLDFGen$Trait," (",QTLDFGen$Threshold, " LOD)")
#Add in dot for sh2 and su1
#su1
#su1 is at chr 4 41,369,510..41,378,299
su1 <- c("Su1",4,round((41369510+41378299)/2), "general",NA,NA,"Su1",NA,NA,NA,NA,NA)
#sh2
#sh2 is at Chr3:216,414,684..216,424,048
sh2 <- c("Sh2",3,(216414684+216424048)/2, "general",NA,NA,"Sh2",NA,NA,NA,NA,NA)
#se1
#se1 is at Chr2:229,983,005..229,983,526
se1 <- c("Se1",2,(229983005+229983526)/2, "general",NA,NA,"Se1",NA,NA,NA,NA,NA)

genes <- data.frame(rbind(rbind(su1,sh2),se1))
colnames(genes) <- colnames(QTLDFGen)
genes$Chrom <- as.numeric(genes$Chrom)
genes$Position <- as.numeric(genes$Position)

#Write those aggregated QTL to file
outfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,DataSet,GWASPolyRunVersion,"_",Thresh,"_SignificantQTL_Aggregated.csv", sep = "")
write.csv(QTLDFGen, outfile)
#get the positions of the snps
genoPosInfo <- geno[,3:4]
#find the final snp position for each chromosome. This is what will determind the size of the chromosome map in the figure 
chrMax <- vector()
for(chr in 1:10){
chrMax[chr] <- max(genoPosInfo$pos[which(genoPosInfo == chr)])}


# make chromosome and end1 as numeric, trait as factor
QTLDFGen$Chrom <- as.numeric(QTLDFGen$Chrom)
QTLDFGen$Position <- as.numeric(QTLDFGen$Position)
QTLDFGen$Trait <- as.factor(QTLDFGen$Trait)

if( length(unique(QTLDFGen$Trait)) != 5){
  clrs <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff")
}
if( length(unique(QTLDFGen$Trait)) == 5){
  clrs <- c("#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff")
}

#create a new table of maize chromosome length and make it a data frame
maize_chromosomes <- cbind(chromosome = c(1:10), start = c(rep(0,10)), end = c(chrMax))
maize_chromosomes <- data.frame(maize_chromosomes)
str(maize_chromosomes)
ytitle <- expression(paste("Genomic positions (",italic(Mb),")"))

# pdf(paste0("Figures/GWASpoly/WSMDP_Carb_GWASpoly_",GWASPolyRunVersion,"_",DataSet,"_AllQTLPosition_GeneralModel.pdf"),width = 10)
tiff(paste0("Figures/GWASpoly/WSMDP_Carb_GWASpoly_",GWASPolyRunVersion,"_",DataSet,"_AllQTLPosition_GeneralModel.tiff"),width = 800)
p <- ggplot(QTLDFGen, aes(Chrom, Position)) +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "black", size = 16, inherit.aes = FALSE) +
  geom_segment(data = maize_chromosomes, aes(x = chromosome, xend = chromosome, y = start, yend = end), lineend = "round", color = "white", size = 15, inherit.aes = FALSE) +
  scale_y_reverse(breaks = seq(3.5e8, 0, -50e6), labels = c(350, seq(300, 0, -50)), limits = c(3.5e8, 0)) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
  ylab(ytitle) + xlab ("Chromosome") +
  theme_classic() +
  theme(axis.text.x = element_text(size=18), axis.title.x = element_text(size= 18), axis.text.y = element_text(size=18), axis.title.y = element_text(size=18)) +
  theme(legend.text = element_text(size=12), legend.title = element_text(size=12)) +
  geom_point(aes(color= Trait), position = position_dodge(width = 0.4), size = 3, alpha = 1) +
  # geom_point(data = genes, aes(x = Chrom, y = Position,color = Trait), size = 5)+data = genes, aes(x = Chrom +.5, y = Position, xend = (Chrom+.1), yend = Position, color = Trait)
  geom_segment(data = genes, aes(x = Chrom +.4, y = Position, xend = (Chrom+.25), yend = Position), size = 1.5,
                 arrow = arrow(length = unit(.25, "cm")))+
  geom_text(data = genes,label= genes$Trait, nudge_x = 0.5, nudge_y = 10000000)+
  # theme(legend.title = element_blank(), legend.position = c(0.50,0.08), legend.direction = "horizontal", legend.text=element_text(size=15)) +  guides(colour = guide_legend(nrow = 2)) +
  scale_color_manual(values = clrs)
  print(p)
  # add_h_arrow(p, genes$Chrom, genes$Position, label = genes$Trait)
dev.off()
return (p)
}

# p <- FullGWASVisualize("EndoFixedEffect_FDRThresh_2021-09-27")
# p <- FullGWASVisualize("EndoFixedEffect_FDRThresh_2021-10-20")
p <- FullGWASVisualize("EndoFixedEffect_FDRThresh_2021-12-09")
p <- p +  theme(legend.position = "none")
# g <- FullGWASVisualize("NoFixedEffect_FDRThresh_2021-10-20")
g <- FullGWASVisualize("NoFixedEffect_FDRThresh_2021-12-09")
g <- g +  theme(legend.title = element_blank(), legend.position = c(0.75,0.2), legend.text=element_text(size=15)) +  guides(colour = guide_legend(nrow = 8)) 
  
# pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoandEndoFixedEffect_FDRThresh_2021-10-20_NFBLUP_AllQTLPosition_GeneralModel.pdf",width = 10, height = 14,)
pdf("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoandEndoFixedEffect_FDRThresh_2021-12-09_NFBLUP_AllQTLPosition_AdditiveModel.pdf",width = 10, height = 14,)
# tiff("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoandEndoFixedEffect_FDRThresh_2021-10-20_NFBLUP_AllQTLPosition_GeneralModel.tiff",width = 1000, height = 1400)
# tiff("Figures/GWASpoly/WSMDP_Carb_GWASpoly_NoandEndoFixedEffect_FDRThresh_2021-12-09_NFBLUP_AllQTLPosition_AdditiveModel.tiff",width = 1000, height = 1400)
figure <- ggarrange(g,p,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2,
                    font.label = list(size = 24))

print(figure)
dev.off()
# FullGWASVisualize("NoFixedEffect_FDRThresh_2021-09-10","WFBLUP")

#########################
###Compare BLUPS by Endosperm with Tukey###
#########################

endoTukey <- function(BlupDF, DataSet){
#establish a file to put the results in
EndoCompareFile <- paste0("Data/OutputtedData/WSMDP_TukeyHSD_CarbBLUP_ComparedbyEndo_",DataSet,"Verbose.txt")
EndoCompareTableFile <- paste0("Data/OutputtedData/WSMDP_TukeyHSD_CarbBLUP_ComparedbyEndo_",DataSet,".csv")
EndoNum <- length(unique(BlupDF$endo))
if( EndoNum == 4){
  LSMHolder <- data.frame("Trait" = blups,  "sh2" = rep("sh2",7), "su1" = rep("su1",7), "su1se1" = rep("su1se1",7), "su1sh2-i" = rep("su1sh2-i",7) )
}
if( EndoNum == 6){
  LSMHolder <- data.frame("Trait" = blups, "sh2" = rep(NA,7), "su1" = rep(NA,7), "su1se1" = rep(NA,7), "su1sh2-i" = rep(NA,7),  "aeduxw" = rep(NA,7), "field" = rep(NA,7) )
}

#iterate through all the traits
for(i in 1:length(blups)){
  #paste in a header
  blup <- blups[i]
  cat(paste0("Compare the ",blup," by endosperm mutant groups."), file=EndoCompareFile, sep="\n", append=TRUE)
  EndoCompare<-aov(BlupDF[,blup]~ BlupDF$endo)
  # EndoCompare<-aov(CleanedInfoNF$Starch~CleanedInfoNF$endo)
  
  
  #find the least square means by group
  LSD<-lsmeans(EndoCompare, ~endo)
  #get the significant differences between the least square means
  LSM<-cld(LSD,Letters = LETTERS, decreasing=T)
  LSMHolder[i,2] <- paste(round(LSM[which(LSM$endo=="sh2"),2],2),"+/-", round(LSM[which(LSM$endo=="sh2"),3],2), LSM[which(LSM$endo=="sh2"),7])
  LSMHolder[i,3] <- paste(round(LSM[which(LSM$endo=="su1"),2],2),"+/-", round(LSM[which(LSM$endo=="su1"),3],2), LSM[which(LSM$endo=="su1"),7])
  LSMHolder[i,4] <- paste(round(LSM[which(LSM$endo=="su1se1"),2],2),"+/-", round(LSM[which(LSM$endo=="su1se1"),3],2), LSM[which(LSM$endo=="su1se1"),7])
  LSMHolder[i,5] <- paste(round(LSM[which(LSM$endo=="su1sh2-i"),2],2),"+/-", round(LSM[which(LSM$endo=="su1sh2-i"),3],2), LSM[which(LSM$endo=="su1sh2-i"),7])
  if(EndoNum == 6){
    LSMHolder[i,6] <- paste(round(LSM[which(LSM$endo=="aeduwx"),2],2),"+/-", round(LSM[which(LSM$endo=="aeduwx"),3],2), LSM[which(LSM$endo=="aeduwx"),7])
    LSMHolder[i,7] <- paste(round(LSM[which(LSM$endo=="field"),2],2),"+/-", round(LSM[which(LSM$endo=="field"),3],2), LSM[which(LSM$endo=="field"),7])
    
  }

  #output that data
  out <- capture.output(LSM)
  cat(out, file=EndoCompareFile, sep="\n", append=TRUE)
}
write.table(LSMHolder, EndoCompareTableFile, sep = ",", row.names = FALSE)
}
endoTukey(NFBlupsGenoJustPheno,"NFBlup")
endoTukey(WFBlupsGenoJustPheno,"WFBlup")

# 
enviTukey <- function(BlupDF, DataSet){
  #establish a file to put the results in
  EnviCompareFile <- paste0("Data/OutputtedData/WSMDP_TukeyHSD_CarbBLUP_ComparedbyEnvi_",DataSet,"Verbose.txt")
  EnviCompareTableFile <- paste0("Data/OutputtedData/WSMDP_TukeyHSD_CarbBLUP_ComparedbyEnvi_",DataSet,".csv")
  EnviNum <- length(unique(BlupDF$Envi))
  if( EnviNum == 4){
    LSMHolder <- data.frame("Trait" = blups,  "2014NewYork" = rep("sh2",7), "2015NewYork" = rep("su1",7), "2014Wisconsin" = rep("su1se1",7), "2015Wisconsin" = rep("su1sh2-i",7) )
  }
  blups <- colnames(BlupDF)[5:11]

  #iterate through all the traits
  for(i in 1:length(blups)){
    #paste in a header
    blup <- blups[i]
    cat(paste0("Compare the ",blup," by Envi groups."), file=EnviCompareFile, sep="\n", append=TRUE)
    #TO FIX this isn't working, frustratingly.
    EnviCompare<-aov(pull(BlupDF[blup])~ BlupDF$Envi)
    # EnviCompare<-aov(CleanedInfoNF$Starch~CleanedInfoNF$Envi)


    #find the least square means by group
    LSD<-lsmeans(EnviCompare, ~Envi)
    #get the significant differences between the least square means
    LSM<-cld(LSD,Letters = LETTERS, decreasing=T)
    LSMHolder[i,2] <- paste(round(LSM[which(LSM$Envi=="2014NewYork"),2],2),"+/-", round(LSM[which(LSM$Envi=="2014NewYork"),3],2), LSM[which(LSM$Envi=="2014NewYork"),7])
    LSMHolder[i,3] <- paste(round(LSM[which(LSM$Envi=="2015NewYork"),2],2),"+/-", round(LSM[which(LSM$Envi=="2015NewYork"),3],2), LSM[which(LSM$Envi=="2015NewYork"),7])
    LSMHolder[i,4] <- paste(round(LSM[which(LSM$Envi=="2014Wisconsin"),2],2),"+/-", round(LSM[which(LSM$Envi=="2014Wisconsin"),3],2), LSM[which(LSM$Envi=="2014Wisconsin"),7])
    LSMHolder[i,5] <- paste(round(LSM[which(LSM$Envi=="2015Wisconsin"),2],2),"+/-", round(LSM[which(LSM$Envi=="2015Wisconsin"),3],2), LSM[which(LSM$Envi=="2015Wisconsin"),7])


    #output that data
    out <- capture.output(LSM)
    cat(out, file=EnviCompareFile, sep="\n", append=TRUE)
  }
  write.table(LSMHolder, EnviCompareTableFile, sep = ",", row.names = FALSE)
}
enviTukey(CleanedInfoNF,"NFCleaned")



justthebitsNF <- NFBlupsGenoJustPheno[2:8]
colnames(justthebitsNF) <- c("Starch","Total poly","WSP","Glucose","Fructose","Sucrose","Total sugar")
pdf(paste("Figures/WSMDP_AllNIRPred_MixedEqn_CorrelationfromPSYCH_BLUPs_NoField.pdf",sep=""))
# tiff(paste("Figures/WSMDP_AllNIRPred_MixedEqn_CorrelationfromPSYCH_BLUPs_NoField.tiff",sep=""))
pairs.panels(justthebitsNF, stars = TRUE, cex = 1,cex.lab = 2, cex.axis = 1)
# cor.plot(justthebitsNF, stars = TRUE, cex.cor = 3.5, cex.axis = 1)
# pairs(justthebitsNF, stars = TRUE, cex.cor = 3.5, cex.axis = 1)
dev.off()
pdf(paste("Figures/WSMDP_EndoInbred_Barchart.pdf",sep=""))
ggplot(WFBlupsGenoJustPheno, aes(x = endo, fill = endo))+
         geom_bar(stat="count")
dev.off()

grep(blups, ".*")

df <- data.frame(Carb = c("Starch", "Total.Polysaccharides", "WSP","Glucose","Fructose","Sucrose", "Total.Sugar"),Heritability =c(0.86,0.94,0.94,0.88,0.86,0.96,0.95) )
pdf(paste("Figures/WSMDP_HeritabilityBarchart.pdf",sep=""),width = 10,height = 10)
ggplot(df, aes(x = Carb,y = Heritability))+
         geom_bar(stat = "identity",fill="steelblue") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
#########################
### LD Visualization ###
#########################


#Import the name of the file. This matrix was calculated on Tassel using the Linkage Disequilibrium tool and the output was saved as a text file
#since the resulting text file was 8GB, it was pruned using AWK on linux to only include the columns for distance and r2 and only those with results, non with NA as the r2
LDMatrixfileName <- "WSMDP_SeqG_NoLD_InbredsPruned_LDDecay_200Window_JustLociDistR2_NoNAs"
fullpath <- paste("Data/RawData/",LDMatrixfileName,".txt",sep = "")
#Read in the .txt LD matrix file
readin <- fread(fullpath, sep='\t')

#Uncomment if all the chromosomes want to be visualized individually, or just the whole geno.
# chrmiterant <- c("All",6)
chrmiterant <- c("All")
for(j in 1:length(chrmiterant)){
  chrmnum = chrmiterant[j]
  data <- readin
  
  #Look at only one chromosome at a time
  if(chrmnum != "All"){
    data <- data[data$Locus1 == chrmnum]}
  
  #uncomment these if the Matrix file has more than two columns. Eliminate any columns not r2 or dist
  # data <- data[,c(13,14)]
  # data <- data[,c(3,4)]
  data <- data[,c(5,6)]
  
  #seperate the distances between chromosomes into 20 bins
  numbins = 20
  #change this if wanting to not log transfrom the data
  log_transform = TRUE
  absolute_distance = TRUE
  
  #change this if you want to eliminate any comparisons less than 500bp apart
  removefirst500 = FALSE
  names(data)=c("dist","rsq")
  data$dist <-as.numeric(data$dist)
  medians=get.median(data, log_transform=log_transform, numbins = numbins, removefirst500 = removefirst500)
  write.table(medians, file="medians.txt",quote=FALSE,sep="\t", col.names=TRUE, row.names=T)
  
  probs=c(0.5, 0.6, 0.7, 0.8, 0.9)
  percents=get.percentiles(data, probs=probs, log_transform=log_transform, numbins = numbins, removefirst500 = removefirst500)
  write.table(percents, file="percents.txt",quote=FALSE,sep="\t", col.names=TRUE, row.names=T)
  mean=tapply(X=data$rsq, INDEX=cut(data$dist, breaks=numbins), FUN=mean, na.rm=T)
  
  percents=read.table('percents.txt',sep='\t',header=T,row.names=1)
  medians=read.table('medians.txt',sep='\t',header=T,row.names=1)
  probs=c('50%','60%','70%','80%','90%')
  colnames(percents) <- probs
  ## Add all together
  row.names(percents)[20]
  
  #the rownames are the starting and end distance of each bin. convert that into a more usable format
  xvals = as.numeric(sub(rownames(percents), pattern="\\((.+),.+", repl="\\1"))
  xvals = c(xvals,8.49)
  mean.xvals <- NULL
  
  #find the middle between each bin border. This will be used for the Xlabel
  for (i in 1:c(length(xvals)-1)){
    mean.xvals <- c(mean.xvals,mean(c(xvals[i],xvals[i+1])))
  }
  
  #determine the file name to save the visualization
  if(removefirst500 == TRUE){
    outfile = paste("Figures/LD_20bins_less500bprmv_from",LDMatrixfileName,"_Chrm",chrmnum,".pdf",sep = "")}
  if(removefirst500 == FALSE){
    outfile = paste("Figures/LD_20bins_from",LDMatrixfileName,"_Chrm",chrmnum,".pdf",sep = "")
  }
  
  #plot and save the visualization
  {pdf(outfile, pointsize=8,family='serif',width = 5,height = 3.5)
  # {jpeg(outfile, pointsize=8,family='serif',width = 500,height = 300)
    lwd=1
    ylim=c(0,1)
    xlim= c(0, range(mean.xvals)[2] * 1.1)	# Increase the upper xlim by a little bit
    par(cex=1,mar=c(5,5,1,1))
    if(log_transform == TRUE){
      xlabel = "Physical distance (log scale)"
    }else{
      xlabel = "Distance"
    }
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlabel, ylab=bquote("Linkage disequilibrium"~(italic(r)^2)), xaxt="n")	#CHD added
    if(log_transform == TRUE){
      x_ticks = 0:8
      x_labels=c("1 bp","10 bp","0.1 kb","1 kb","10 kb","100 kb","1 Mb","10 Mb","100 Mb") #Jason original
    } else{
      x_ticks = seq(from=xlim[1], to=xlim[2], length.out=5)
      x_labels = x_ticks
    }
    axis(side=1, at=x_ticks, labels=x_labels)
    for(j in 1:length(probs)){
      lines(x=mean.xvals, y=percents[,j], lwd=lwd, col='black', lty=j)
    }
    text(x=sort(xvals)[8], y=c(0.11,0.30,0.46,0.68,0.9), labels=colnames(percents), pos=4, cex = 1)
    dev.off()
  }
}


#########################
### Candidate Gene Discovery ###
#########################

## Identify candidate genes using GenomicRanges

# Prep the genes ----------------------------------------------------------
library(rtracklayer)
gff = import(gzfile('Data/RawData/ZmB73_5b.60_FGS.gff3.gz'))
# gff = import(gzfile('Data/RawData/Zea_mays.AGPv3.22.gff3.gz'))
head(gff)


# Create tibble from gff file removing non-chromosomal mitochondrial, and chloroplast genes
library(data.table)
library(tidyverse)

gffFiltered = as_tibble(gff) %>% 
  # select(seqnames, start, end, type, ID, description) %>% 
  filter(type == "gene") %>% 
  # filter(!grepl("B73",seqnames)) %>% 
  # filter(grepl("Chr*",seqnames)) %>% 
  # filter(!grepl("Mt",seqnames)) %>%
  filter(!grepl("chrMt",seqnames)) %>%
  # filter(!grepl("Cp",seqnames)) %>%
  # filter(!grepl("UNMAPPED",seqnames))%>%
  filter(!grepl("chrUNKNOWN",seqnames))%>%
  # filter(!grepl("scaffold",seqnames)) %>%
  filter(!grepl("chrPt",seqnames)) %>%
  separate(seqnames, into = c("Intro","seqnames"), sep = "r") %>%
  mutate(seqnames = as.numeric(seqnames)) %>%
  arrange(seqnames, start)

head(gffFiltered)



## Gene ranges
library(GenomicRanges)
gR = with(gffFiltered, GRanges(seqnames = seqnames,
                       ranges = IRanges(start = start, end = end),
                       ids = ID,
                       description = biotype))
                       # description = description))

## Read in SNP results should have SNP, CHROM, POS, p.value, Trait columns
n_marks = 8e4

#Pull from this file
GWASPolyRunVersion <- "EndoFixedEffect_FDRThresh_2021-09-27"
# GWASPolyRunVersion <- "NoFixedEffect_FDRThresh_2021-09-27"
DataSet = "NFBLUP"
Thresh = "FDR"
outfile <- paste0("WSMDP_Carb_GWASpoly_",Seq,DataSet,GWASPolyRunVersion,"_",Thresh,"_SignificantQTL_Aggregated")
GWASResultsFile <- paste("Data/OutputtedData/GWASpoly/",outfile,".csv", sep = "")
snpDat = fread(GWASResultsFile, data.table = FALSE)
#ToOnlyLookAtSigLoci
# snpDatFull = fread(GWASResultsFile, data.table = FALSE)
# snpDat <- snpDatFull[which(snpDatFull$pval < 0.05),]

head(snpDat)
bonf = 0.05
window = 250e3

## SNP ranges
sR = snpDat %>% 
  # filter(Method == "FarmCPU") %>% 
  # filter(pval<=bonf) %>%
  arrange(Chrom ,Position) %>%
  with(., GRanges(seqnames = Chrom, 
                  ranges = IRanges(start = Position - window, end = Position + window),
                  pheno = Trait, 
                  SNP = Marker))

## Find overlaps
#This somehow finds things of same positions regardless of chrom number, so thats no good
sR@seqnames
gR@seqnames
my_genes = findOverlaps(gR, sR)
# split(gR, seqnames(gR))$'1'
# 
# 
# my_genes = findOverlaps(gR@ranges, sR@ranges)
# my_genes = as_tibble(findOverlaps(split(gR, seqnames(gR))$'1'@ranges, split(sR, seqnames(sR))$'1'@ranges))
# my_genes <- my_genes %>% add_column("Chrm" = 1)
# 
# for(i in 2:10)
# {
#   tempgR <- split(gR, seqnames(gR))[i]
#   tempsR <-split(sR, seqnames(sR))[i]
#   tempGenes <- as_tibble(findOverlaps(tempgR@unlistData@ranges, tempsR@unlistData@ranges))
#   tempGenes <- tempGenes %>% add_column("Chrm" = i)
#   my_genes = rbind(my_genes, tempGenes)
# }

## Index variables
gmind <- my_genes@from
smind <- my_genes@to

#this is donked
## Create a table from genomic ranges and SNP ranges
genes_results = tibble(GeneID = gR@elementMetadata@listData$ids[c(gmind)],
                       SNP = c(sR@elementMetadata@listData$SNP[smind]),
                       Phenotype = c(sR@elementMetadata@listData$pheno[smind])) %>%
  inner_join(., gffFiltered, by = c("GeneID" = "ID")) %>%
  arrange(Phenotype, seqnames, start)
write.table(genes_results,paste0("Data/OutputtedData/GWASpoly/",outfile,"_CandidateGenes_B73v2.csv"))
