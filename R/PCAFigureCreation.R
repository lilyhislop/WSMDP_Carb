
PCAFigureCreation <- function(PCA,pc.perc,genoinfo,infilename, Colourant ){
  #ForDebugging
  genoinfo <- matching
  Colourant <- "Program"
progfile <- paste("Figures/Plot_PCA_",infilename,"_by_",Colourant,".png",sep = "")

colors <-c("#f0e442","#000000", "#cc79a7","#0072b2", "#009e73","#d55e00")
# symbols <- c(0,1,15,2,3,4,18,6)
symbols <- c(25,16,17,18,19,15)

if(Colourant == "Program"){
  #label with program
  genoinfo$ProgramOthered <- genoinfo$Program
  toOther <- c("FL","GA","MA","ME","MO","ND","NJ","OH","SC","SP","TN","USDA","")
  levels(genoinfo$ProgramOthered) <- c(levels(genoinfo$ProgramOthered),"Other")
  genoinfo$ProgramOthered[genoinfo$Program %in% toOther]<- "Other"
  colors <- c("slategray3", "violet", "royalblue4", "palegreen", "mediumpurple", "khaki", "grey", "red")
  symbols <- c(0,1,15,2,3,4,6,18)
  Colourant <- "ProgramOthered"
}


if(Colourant == "SusceptibilityRating03"){
  colors <- c("royalblue4", "palegreen", "mediumpurple", "red")
symbols <- c(0,1,15,2)}

png(progfile,width = 750,height = 750)
#label PCA with region
tab <- data.frame(sample.id = PCA$sample.id, 
                  pop = factor(genoinfo[,Colourant])[match(genoinfo$GenoName,PCA$sample.id)],
                  EV1 = PCA$eigenvect[,1],
                  EV2 = PCA$eigenvect[,2],
                  stringsAsFactors = FALSE)
plot(tab$EV2, tab$EV1,pch = symbols[as.integer(tab$pop)], col = colors[as.integer(tab$pop)], xlab = paste("PC 2 (",(round(pc.perc,1))[2],"%)" ,sep = ""), ylab = paste("PC 1 (",(round(pc.perc,1))[1],"%)" ,sep = ""))
legend("bottomleft", legend = levels(tab$pop), pch = symbols, col = colors)
dev.off()
}



