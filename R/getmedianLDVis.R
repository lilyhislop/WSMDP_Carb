#this is a funtion to calculate what's the median LD the compared SNPS
get.median=function(x, probs=5, log_transform=log_transform, numbins=numbins, removefirst500 = removefirst500){
  mydists = x$dist
  if(removefirst500 == TRUE){
    x <- x[-which(mydists < 501),]
    mydists <- mydists[-which(mydists < 501)]}
  if(log_transform == TRUE){
    is_negative = mydists<0 & !is.na(mydists) 
    mydists = abs(mydists)	# Have to abs-value b/c log doesn"t do negatives
    print(range(mydists))
    mydists[mydists==0]=1 # To avoid problems with log-transforming 0-values
    mydists = log10(mydists)	
    mydists[is_negative] <- mydists[is_negative] * -1	# Restore to negative to indicate is before the gene
  }
  results=as.list(tapply(X=x$rsq, INDEX=cut(mydists, breaks=numbins), FUN=median, na.rm=T))
  results=do.call(what=rbind, args=results)
  return(results)
}