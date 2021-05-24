
##########Linear Model Analysis!##########
library(lme4)
library(lattice)

#establish working directory
setwd("C:/Users/LHislop/Documents/GitHub/WSMDP_Carb")
CleanedInfoWF <- read.csv(file = "Data/OutputtedData/CleanedInfoWFOutput.csv", stringsAsFactors  = TRUE )

#A differently callibrated version of the same dataset
# CleanedInfoNF <- read.csv(file = "Data/OutputtedData/CleanedInfoNFOutput.csv")


head(CleanedInfoWF)
colnames(CleanedInfoWF)

########## lm version ##########


#establish dataframe to store variances from each factor
VarDF <- data.frame("Carb" = colnames(CleanedInfoWF)[6:12],"Geno" = rep(NA,7),"Envi" = rep(NA,7), "Superblock" = rep(NA,7),"Block" = rep(NA,7),"Geno:Envi" = rep(NA,7),"Residuals"= rep(NA,7))

#for each carb version, look what factors are incluencing the variation
for(i in 1:7){
  #formula is carb ~ check+ Gene + Envi + Gene*Envi + superblock(envi) + block(superblockxenvi)+ error 
  formula1 <- paste(colnames(CleanedInfoWF)[i+4],"~ Check + Variety*Envi + superblock:Envi + block:superblock:Envi")
  fit1 <- lm(formula1,data=CleanedInfoWF)
  AIC1 <- extractAIC(fit1)
  AIC1
  summary(fit1)$r.square
  print(anova(fit1))
  out <- anova(fit1)
  SStotal <- sum(out$`Sum Sq`)
  for(j in 1:6){
    #variance explained is caluclated by the sum of squares divided by the sum of squares total
    VarDF[i,j+1] <- (out$`Sum Sq`[j]/SStotal)}
}

VarDFMelt <- melt(VarDF)

#Visualization of percent variance for each factor 
barchart(~value|variable, group = factor(Carb), data= VarDFMelt,main = "Percent Phenotypic Variance Explained",layout = c(6,1),
         key = simpleKey(text = colnames(CleanedInfoWF)[5:11],
                         rectangles = TRUE, points = FALSE, space = "right"))



########## lmer version ##########

carbs <- colnames(CleanedInfoWF)[6:12]
for(carb in carbs){
  normality <- shapiro.test(CleanedInfoWF[,carb])
  out <- capture.output(normality)
  
  
  model <- lmer(get(carb) ~ Check + (1|Variety) + Envi + Variety:Envi + superblock/Envi + block/Envi, data=CleanedInfoWF, REML = TRUE)
  summary(model, correlation = FALSE)
  random_effects <- ranef(model)
}

########## Playing down here ##########
DFSubset <- subset(CleanedInfoWF, !is.na(IsExperimental))
str(DFSubset)
colnames(DFSubset)

DFSubset$superblock <- as.factor(DFSubset$superblock)
DFSubset$Col <- as.factor(DFSubset$Col)
DFSubset$Row <- as.factor(DFSubset$Row)
DFSubset$Year <- as.factor(DFSubset$Year)
DFSubset$Check <- as.factor(DFSubset$Check)
DFSubset$block <- as.factor(DFSubset$block)
DFSubset$Rep <- as.factor(DFSubset$Rep)
DFSubset$PlotNum <- as.factor(DFSubset$PlotNum)

#use likelihood ratio test: to test specific fixed effects
#useAIC to test model
#with block, superblock as stand alone fixed effects, envi as random
model <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + block, data=DFSubset, REML = TRUE)
summary(model)
plot(model)
qqmath(ranef(model, condVar=TRUE))
#looks pretty ok, probably don't need tranformation
qqnorm(residuals(model))
extractAIC(model)
AIC(model)


#with block, superblock, row, column as stand alone fixed effects, envi as random
modelRC <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + Row + Col + block, data=DFSubset, REML = TRUE)
summary(modelRC)
plot(modelRC)
qqmath(ranef(modelRC, condVar=TRUE))
qqnorm(residuals(modelRC))
extractAIC(modelRC)
AIC(modelRC)

#with block, superblock, row,  as stand alone fixed effects, envi as random
modelR <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi) + superblock + Row +  block, data=DFSubset, REML = TRUE)
summary(modelR)
plot(modelR)
qqmath(ranef(modelR, condVar=TRUE))
qqnorm(residuals(modelR))
extractAIC(modelR)
AIC(modelR)


#with block, superblock, row, column  effects and envi all as random
modelRCNest <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi/superblock/block) + (1|Envi:Row) + (1|Envi:Col) , data=DFSubset, REML = TRUE)
summary(modelRCNest)
plot(modelRCNest)
qqmath(ranef(modelRCNest, condVar=TRUE))
qqnorm(residuals(modelRCNest))
extractAIC(modelRCNest)
AIC(modelRCNest)


#with block, superblock nested random effect
modelnest <- lmer(Total.Sugar ~ Check + (1|Variety) + (1|Variety:Envi) + (1|Envi/superblock/block), data=DFSubset, REML = TRUE)
summary(modelnest)
plot(modelnest)
qqmath(ranef(modelnest, condVar=TRUE))
#looks pretty ok, probably don't need tranformation
qqnorm(residuals(modelnest))
extractAIC(modelnest)
AIC(modelnest)



#model workied out with raegan
modelRae<- lmer(Total.Sugar ~  (1|superblock:Variety) + (1+Envi|Variety) + (1|superblock/block), 
                data=DFSubset, REML = TRUE)
summary(modelRae)
plot(modelRae)
qqmath(ranef(modelRae, condVar=TRUE))
qqnorm(residuals(modelRae))
extractAIC(modelRae)
AIC(modelRae)








#if there is no evidence of a difference, go for the simpler model. If P>0.05, then there is no evidence that they are meaningfully different and just go with the simpler one
# p > 0.05, theres no evidence against the null H0 that sigma^2_whatever = 0 so youprefer the simpler model
anova(modelnest,modelRCNest)

anova(model,modelnest)

anova(model,modelRC)

# Questions
# How do you justify it being better than the published model
# whats normal in the field. If not a lot of copies its fixed. Ask reagan about whats normal
# If the published model isn't as good of a fit, but how different are the estimates? Look at the Ranef
# if the ranef isn't changing much, maybe include just because its a narrative that makes sense
