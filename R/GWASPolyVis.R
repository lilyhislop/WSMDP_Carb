GWASPolyVis <- function(GWASPolyRunVersion,trait,data3,Seq,DataSet, Thresh = "FDR"){
  
  visfileprefix <- paste0("Figures/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,"_",DataSet,"_",GWASPolyRunVersion,"_",trait,"_",Thresh)
  
  QQplotfile <- paste(visfileprefix,"_QQplot_General.png", sep = "")
  png(QQplotfile, width = 500, height = 500)
  #,"1-dom-alt-alt","1-dom-alt-ref","1-dom-ref-alt","1-dom-ref-ref")
  print(qq.plot(data3,trait=trait,model="general"))
  
  dev.off()
  
  data4 <- set.threshold(data3,method=Thresh,level=0.05)
  #Here's the QTLS found
  print(get.QTL(data4))
  data5 <- get.QTL(data4) 
  # knitr::kable(data5)
  # Models <- data5$Model
  # Markers <- data5$Marker[which(data5$Score >= data5$Threshold)] ## this one works!!!
  # trait <- data5$Trait
  

  # knitr::kable(fit.ans,digits=3)
  
  # p <- LD.plot(data4)
  # p + xlim(0,30) 
  # 
  # system.time(beep(for (i in 1:nrow(data5)) {
  #   R2 <- fit.QTL(data=data4,trait=trait[i],qtl=data.frame(Q1=c(Markers[i], Models[i])))
  #   data5[i,]$R2 <- R2[,3]
  #   print(paste(trait[i], i, "of", nrow(data5))) # "of 97")
  #   print(data5[i,])
  # }, sound = "fanfare")) # takes 10 minutes to print everything <3 <3 <3
  # 
  # data5

  # print(fit.QTL(data2, trait))
  
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

  MANAddplotfile <- paste(visfileprefix,"_Manhattan_AdditiveModel.png", sep = "")
  png(MANAddplotfile,width = 750, height = 500)
  print({
  par(mfrow=c(1,1))
  manhattan.plot(data4,trait=trait,model="additive")
  })
  dev.off()

  Scoresfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,DataSet,GWASPolyRunVersion,"_",trait,"_",Thresh,"_scores.csv", sep = "")
  Effectfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,DataSet,GWASPolyRunVersion,"_",trait,"_",Thresh,"_effects.csv", sep = "")
  write.GWASpoly(data4, trait, filename=Scoresfile, what = "scores", delim = ",")
  write.GWASpoly(data4, trait, filename=Effectfile, what = "effects", delim = ",")
  data5File <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,DataSet,GWASPolyRunVersion,"_",trait,"_",Thresh,"_SignificantQTL.csv", sep = "")
  write.table(data5,
              append = FALSE,
              file = data5File,
              sep = ",",
              dec = ".",
              row.names = FALSE,
              col.names = TRUE)
  #  if(length(data5)>0){
  #   fit.ans <- fit.QTL(data=data4,trait=trait,
  #                      qtl=data5[,c("Marker","Model")],
  #                      fixed=data.frame(Effect="endo",Type="factor"))
  # QTLfile <- paste("Data/OutputtedData/GWASpoly/WSMDP_Carb_GWASpoly_",Seq,DataSet,GWASPolyRunVersion,"_",trait,"_",Thresh,"_QTLswithEffects.csv", sep = "")
  # write.table(fit.ans,
  #             append = FALSE,
  #             file = QTLfile,
  #             sep = ",",
  #             dec = ".",
  #             row.names = FALSE,
  #             col.names = TRUE)}

}
