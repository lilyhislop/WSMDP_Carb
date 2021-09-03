##########################
### Hislop GWASpoly Error Replication####
### With smaller datasets and simplified code for easy investigation!####
### 2021.09.03 #########
########################


###### So the problem is that I have 2 traits. "Total.Polysaccharides.BLUP" is throwing an error when i run their GWAS Results through fit.QTL. "WSP.BLUP" is not.
# This is an abbreviated replicated version of the problem I am having with my full dataset

#########################
### GWASPoly ###
#########################
#Establish Libraries and working directory
library("devtools")#for intsalling from github
install_github("jendelman/GWASpoly")
library(GWASpoly)#for running Gwas
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

#Read in packages
source("R/GWASPolyVis.R")
source("R/GWASPolyRunner.R")
source("R/hmpToNumeric.R")
source("R/WritePhenoGenotoFile.R")

#read in genetic info post MAF
Seq <- "SeqG"
hmppath <- paste0("Data/WSMDP_",Seq,".hmp.txt")

geno <- fread(hmppath,skip = "rs#")
colnames(geno)<-gsub(colnames(geno), pattern = ":.*", replacement = "")

NFBlups <- read.table("Data/CarbPhenoInformationforDebug.txt")
NFBlups$endo <- as.factor(NFBlups$endo)

c1 <- which(colnames(NFBlups)=="Total.Polysaccharides.BLUP")
c2 <- which(colnames(NFBlups)=="WSP.BLUP")
blups <- colnames(NFBlups[c2:c1])

for(blup in blups){
  GWASPolyRunner(NFBlups,geno,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),Seq,"NFBLUP","endo","factor")
}
