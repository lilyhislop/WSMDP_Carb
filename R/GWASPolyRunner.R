GWASPolyRunner <- function(phenoSubsetGeno,geno,trait, RunName, Seq, DataSet,fixedEffect = NULL,fixedType = NULL, Thresh = "FDR"){
# #DEBugging Info
# phenoSubsetGeno <- NFBlupsGenoJustPheno
# trait <- blup
# RunName <- paste0("EndoFixedEffect_FDRThresh_",Sys.Date())
# DataSet <- "NFBLUP"
# fixedEffect <- "endo"
# fixedType <- "factor"
#  Thresh = "FDR"

GWASPolyRunVersion <- RunName

#read in hmp data, output numeric format that GWASPoly likes
ifelse(!is.numeric(geno[15,30]),SCMVPanel_nwithpos <- hmpToNumeric(geno),SCMVPanel_nwithpos <- geno)
outfiles <- WritePhenoGenoToFile(GWASPolyRunVersion,trait,phenoSubsetGeno,SCMVPanel_nwithpos,Seq,DataSet)

#now we run the GWASpoly with the files in the proper format
data <- read.GWASpoly(ploidy=2, 
                      pheno.file = outfiles[1],
                      geno.file=outfiles[2],
                      format="numeric",
                      n.traits=7,
                      delim=",")

data2 <- set.K(data, LOCO=TRUE)
params <- set.params(fixed=fixedEffect, fixed.type=fixedType,n.PC = 3, MAF = 0.025) #no fixed effects, MAF should do nothing as it's already been filtered

data3 <- GWASpoly(data2,models=c("general", "additive"),traits=trait, params=params)

#visualize the gwas results
GWASPolyVis(GWASPolyRunVersion, trait, data3,Seq,DataSet, Thresh)
}
