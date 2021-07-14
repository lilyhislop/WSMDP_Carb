hmpToNumeric <- function(hmpData){
  library("scrime")#for recodeSNPS
  SCMVPanel_un <- hmpData
  
  #assemble the info about marker, chromosome and position
  SCMVPanel_map <- SCMVPanel_un[,c('rs#','chrom','pos')]
  colnames(SCMVPanel_map) <- c('Marker','Chrom','Position')
  
  #convert from hapmap to numeric
  noMeta <- SCMVPanel_un[,-(1:11)]
  SCMVMatrix <- as.matrix(noMeta)
  SCMVMatrix_clean <- SCMVMatrix
  SCMVMatrix_clean[which(SCMVMatrix == "A")] <- "AA"
  SCMVMatrix_clean[which(SCMVMatrix == "C")] <- "CC"
  SCMVMatrix_clean[which(SCMVMatrix == "G")] <- "GG"
  SCMVMatrix_clean[which(SCMVMatrix == "T")] <- "TT"
  SCMVMatrix_clean[which(SCMVMatrix == "R")] <- "AG"
  SCMVMatrix_clean[which(SCMVMatrix == "Y")] <- "CT"
  SCMVMatrix_clean[which(SCMVMatrix == "S")] <- "CG"
  SCMVMatrix_clean[which(SCMVMatrix == "W")] <- "AT"
  SCMVMatrix_clean[which(SCMVMatrix == "K")] <- "GT"
  SCMVMatrix_clean[which(SCMVMatrix == "M")] <- "AC"
  SCMVMatrix_clean[which(SCMVMatrix == "N")] <- NA
  SCMVMatrix_clean[which(SCMVMatrix == "0")] <- NA
  SCMVMatrix_clean[which(SCMVMatrix == "-")] <- NA
  
  #some of the SNPs have more than 2 homozygous genotypes. Get rid of them 
  toomanyHomos <- which(apply(SCMVMatrix_clean, 1, function(x) length(unique(x)))>4)
  if(length(toomanyHomos) > 0){
  SCMVMatrix_clean <- SCMVMatrix_clean[-toomanyHomos,]
  SCMVPanel_map<- SCMVPanel_map[-toomanyHomos,]}
  
  #convert to numeric
  SCMVPanel_n <- recodeSNPs(SCMVMatrix_clean, first.ref = FALSE, geno = 0:2, snp.in.col = FALSE)
  str(SCMVPanel_n)
  
  
  #bind numeric info with mapping into
  SCMVPanel_nwithpos <- cbind(SCMVPanel_map, SCMVPanel_n)
  return(SCMVPanel_nwithpos)
}
