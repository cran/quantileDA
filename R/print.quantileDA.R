print.quantileDA<-function(x,...)
{
	out=x
	message("")
	message("Results of the quantile classifier")
	message("")
	cat(paste("Minimum misclassification rates attained at theta = ",paste(round(out$theta.choice,2),collapse = " "),sep=""),"\n")
	cat(paste("Miclassification rate in the training set: ", round(out$me.train,2),sep=""),"\n")
	if (!is.null(out$me.test)) cat(paste("Miclassification rate in the test set: ",round(out$me.test,2),sep=""),"\n")
	
}