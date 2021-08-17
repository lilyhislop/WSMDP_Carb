WritePhenoGenoToFile <- function(GWASPolyRunVersion,trait, phenoSubsetGeno,SCMVPanel_nwithpos,filename,adendum){
  
  #make the necessary phenotype/genotype files, required by the GWASPoly manual
  phenoSubsetGeno_fullpheno <- phenoSubsetGeno
  phenoSubsetGeno_trait_Only <- phenoSubsetGeno[,c('GenoName',trait)]
  
#read out the phenotype info into a format that GWASpoly will like
phenooutfilefull <- paste("Data/WSMDP_Carb_GWASpoly_FullPheno_",GWASPolyRunVersion,Sys.Date(),".csv",sep = "")
write.table(phenoSubsetGeno_fullpheno,
            append = FALSE,
            file = phenooutfilefull,
            sep = ",",
            dec = ".",
            row.names = FALSE,
            col.names = TRUE)

#read out the phenotype info into a format that GWASpoly will like. This time with only susceptibility info incase the other info is confusing
phenooutfiletrait <- paste("Data/WSMDP_Carb_GWASpoly_Pheno_only_",trait,"_",GWASPolyRunVersion, Sys.Date(),".csv",sep = "")
write.table(phenoSubsetGeno_trait_Only,
            append = FALSE,
            file = phenooutfiletrait,
            sep = ",",
            dec = ".",
            row.names = FALSE,
            col.names = TRUE)

#export as the geno info as CSV
# genooutput <- paste("Data/RawData/",filename,adendum,"_",GWASPolyRunVersion,"_numericFormat",Sys.Date(),".csv",sep = "")
genooutput <- paste("Data/RawData/",filename,"_numericFormat.csv",sep = "")
write.table(SCMVPanel_nwithpos,
            append = FALSE,
            file = genooutput,
            sep = ",",
            dec = ".",
            row.names = FALSE,
            col.names = TRUE)
return(c(phenooutfilefull,genooutput))
}
