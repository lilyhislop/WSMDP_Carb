#this is a funtion to calculate what proportion of the compared snps at a given distance apart are in what LD
get.percentiles=function(x, probs=0.95, log_transform=log_transform, numbins=numbins, removefirst500 = removefirst500){
  # 	cat("Numbins is",numbins,"and range of distances is",range(x$dist),"\n")
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
    mydists[is_negative] = mydists[is_negative] * -1	# Restore to negative to indicate is before the gene
  }
  results=as.list(tapply(X=x$rsq, INDEX=cut(mydists, breaks=numbins), FUN=quantile, na.rm=T, probs=probs))
  results=do.call(what=rbind, args=results)
  return(results)
}

