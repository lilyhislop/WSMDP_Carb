GWASPolyVis <- function(GWASPolyRunVersion,trait,data3,filename,adendum){
  visfileprefix <- paste0("Figures/GWASpoly/WSMDP_Carb_GWASpoly_",filename,"_",adendum,"_",GWASPolyRunVersion,"_",trait)
  #look at the QQplots for each tested model
  # QQplotfile <- paste(visfileprefix,"_QQplot.png", sep = "")
  # png(QQplotfile)
  # print({
  #   par(mfrow=c(2,1))
  # models <- c("additive","general")#,"1-dom-alt-alt","1-dom-alt-ref","1-dom-ref-alt","1-dom-ref-ref")
  # #ONE TRAIT AT A TIME
  # for (i in 1:2) {
  #   qq.plot(data3,trait=trait,model=models[i])
  # }
  # 
  #    })
  # dev.off()
  # 
  # QQplotfile <- paste(visfileprefix,"_QQplot_Additive.png", sep = "")
  # png(QQplotfile)
  # print({
  # par(mfrow=c(1,1))
  # models <- c("additive","general")#,"1-dom-alt-alt","1-dom-alt-ref","1-dom-ref-alt","1-dom-ref-ref")
  # qq.plot(data3,trait=trait,model=models[1])
  # 
  # # })
  # dev.off()
  
  QQplotfile <- paste(visfileprefix,"_QQplot_General.png", sep = "")
  png(QQplotfile)
  
  #,"1-dom-alt-alt","1-dom-alt-ref","1-dom-ref-alt","1-dom-ref-ref")
  qq.plot(data3,trait=trait,model="general")
  dev.off()
  
  data4 <- set.threshold(data3,method="permute",level=0.05)
  #Here's the QTLS found
  print(get.QTL(data4))
  
  #Heres the manhattan plots for this
  # MANplotfile <- paste(visfileprefix,"_Manhattan_AllModels.png", sep = "")
  # png(MANplotfile,width = 750, height = 500)
  # par(mfrow=c(2,1))
  # print({
  # manhattan.plot(data4,trait=trait,model="additive")
  # manhattan.plot(data4,trait=trait,model="general")
  # # manhattan.plot(data4,trait=trait,model="1-dom-alt-ref")
  # # manhattan.plot(data4,trait=trait,model="1-dom-alt-alt")
  # # manhattan.plot(data4,trait=trait,model="1-dom-ref-alt")
  # # manhattan.plot(data4,trait=trait,model="1-dom-ref-ref")
  # 
  # })
  # dev.off()
  # 
  
  MANGenplotfile <- paste(visfileprefix,"_Manhattan_GeneralModel.png", sep = "")
  png(MANGenplotfile,width = 750, height = 500)
  print({
  par(mfrow=c(1,1))
  manhattan.plot(data4,trait=trait,model="general")
  # abline(a = get.QTL(data4, model="general")$Threshold[1], b = 0, lty = "dashed")
  })
  dev.off()
  # 
  # MANAddplotfile <- paste(visfileprefix,"_Manhattan_AdditiveModel.png", sep = "")
  # png(MANAddplotfile,width = 750, height = 500)
  # print({
  # par(mfrow=c(1,1))
  # manhattan.plot(data4,trait=trait,model="additive")
  # })
  # dev.off()
  
  Scoresfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",filename,adendum,GWASPolyRunVersion,"_",trait,"_scores.csv", sep = "")
  Effectfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",filename,adendum,GWASPolyRunVersion,"_",trait,"_effects.csv", sep = "")
  write.GWASpoly(data4, trait, filename=Scoresfile, what = "scores", delim = ",")
  write.GWASpoly(data4, trait, filename=Effectfile, what = "effects", delim = ",")
  
  QTLfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",filename,adendum,GWASPolyRunVersion,"_",trait,"_QTLs.csv", sep = "")
  write.table(get.QTL(data4),
              append = FALSE,
              file = QTLfile,
              sep = ",",
              dec = ".",
              row.names = FALSE,
              col.names = TRUE)
}
