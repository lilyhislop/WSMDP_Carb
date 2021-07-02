#########################
#This code was written by Lillian Hislop
#2021.06.23
#WSMDP Carbohydrate Analysis 
#GWAS Degbugging
#########################
library("devtools")
install_github("jendelman/GWASpoly")
library(GWASpoly)#for running Gwas

getwd()
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")

#The hypothetical functions
# GWASPolyRunner(WFBlupsGenoJustPheno[,1:8],geno_scmv,blup,paste0("NoFixedEffect_FDRThresh_",Sys.Date()),"SeqB","WFBLUP")
# GWASPolyRunner(WFBlupsGenoJustPheno,geno_scmv,blup,paste0("EndoFixedEffect_FDRThresh_",Sys.Date()),"SeqB","WFBLUP","endo","factor")


#From GWASPolyRunner Function
filename <- "SeqB"
adendum <- "WFBLUP"
trait <- "Starch.BLUP"

#No Fixed effect
# GWASPolyRunVersion <- "NoFixedEffect_FDRThresh_2021-06-23"
# fixedeffect <- NULL
# fixedType <- NULL

# #With EndoMutant as Fixed Effect
GWASPolyRunVersion <- "EndoFixedEffect_FDRThresh_2021-06-23"
fixedEffect <- "endo"
fixedType <- "factor"

outfiles1 <- paste("Data/WSMDP_Carb_GWASpoly_FullPheno_",GWASPolyRunVersion ,".csv",sep = "")
outfiles2 <- paste("Data/RawData/",filename,adendum,"_",GWASPolyRunVersion,"_numericFormat.csv",sep = "")

#now we run the GWASpoly with the files in the proper format
readGPOut <- read.GWASpoly(ploidy=2, 
                      pheno.file = outfiles1,
                      geno.file=outfiles2,
                      format="numeric",
                      n.traits=7,
                      delim=",")

setKOut <- set.K(readGPOut, LOCO=FALSE)

paramsOut <- set.params(geno.freq = 1 - 10/957,fixed=fixedEffect, fixed.type=fixedType, n.PC = 3, MAF = 0.005) #no fixed effects, MAF should do nothing as it's already been filtered

# THIS IS WHERE THE PROBLEM IS IF USING FIXED EFFECTS. Error: Error in cbind(X, model.matrix(~x, data.frame(x = xx))[, -1]) : 
#   number of rows of matrices must match (see arg 2)
GPOut <- GWASpoly(setKOut,models="general",traits=trait, params=paramsOut)


#From My GWASPolyVis function
visfileprefix <- paste0("Figures/GWASpoly/WSMDP_Carb_GWASpoly_",filename,"_",adendum,"_",GWASPolyRunVersion,"_",trait)

QQplotfile <- paste(visfileprefix,"_QQplot_General.png", sep = "")
png(QQplotfile)
#,"1-dom-alt-alt","1-dom-alt-ref","1-dom-ref-alt","1-dom-ref-ref")
qq.plot(GPOut,trait=trait,model="general")
dev.off()

threshOut <- set.threshold(GPOut,method="FDR",level=0.05)

#Here's the QTLS found
print(get.QTL(threshOut))



MANGenplotfile <- paste(visfileprefix,"_Manhattan_GeneralModel.png", sep = "")
png(MANGenplotfile,width = 750, height = 500)
print({
  par(mfrow=c(1,1))
  manhattan.plot(threshOut,trait=trait,model="general")
  # abline(a = get.QTL(threshOut, model="general")$Threshold[1], b = 0, lty = "dashed")
})
dev.off()
# 
# MANAddplotfile <- paste(visfileprefix,"_Manhattan_AdditiveModel.png", sep = "")
# png(MANAddplotfile,width = 750, height = 500)
# print({
# par(mfrow=c(1,1))
# manhattan.plot(threshOut,trait=trait,model="additive")
# })
# dev.off()

Scoresfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",filename,adendum,GWASPolyRunVersion,"_",trait,"_scores.csv", sep = "")
Effectfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",filename,adendum,GWASPolyRunVersion,"_",trait,"_effects.csv", sep = "")
write.GWASpoly(threshOut, trait, filename=Scoresfile, what = "scores", delim = ",")
write.GWASpoly(threshOut, trait, filename=Effectfile, what = "effects", delim = ",")

QTLfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",filename,adendum,GWASPolyRunVersion,"_",trait,"_QTLs.csv", sep = "")
write.table(get.QTL(threshOut),
            append = FALSE,
            file = QTLfile,
            sep = ",",
            dec = ".",
            row.names = FALSE,
            col.names = TRUE)

