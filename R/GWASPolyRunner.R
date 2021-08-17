GWASPolyRunner <- function(phenoSubsetGeno,geno,trait,RunName,filename,adendum,fixedeffect = NULL,fixedType = NULL){
  GWASPolyRunVersion <- RunName

#read in hmp data, output numeric format that GWASPoly likes
ifelse(!is.numeric(geno[15,30]),SCMVPanel_nwithpos <- hmpToNumeric(geno),SCMVPanel_nwithpos <- geno)
outfiles <- WritePhenoGenoToFile(GWASPolyRunVersion,trait,phenoSubsetGeno,SCMVPanel_nwithpos,filename,adendum)

#now we run the GWASpoly with the files in the proper format
data <- read.GWASpoly(ploidy=2, 
                      pheno.file = outfiles[1],
                      geno.file=outfiles[2],
                      format="numeric",
                      n.traits=7,
                      delim=",")

data2 <- set.K(data, LOCO=TRUE)
params <- set.params(fixed=fixedeffect, fixed.type=fixedType,n.PC = 3, MAF = 0.005) #no fixed effects, MAF should do nothing as it's already been filtered

data3 <- GWASpoly(data2,models=c("diplo-general", "diplo-additive"),traits=trait, params=params)

#visualize the gwas results
GWASPolyVis(GWASPolyRunVersion, trait, data3,filename,adendum)
}
