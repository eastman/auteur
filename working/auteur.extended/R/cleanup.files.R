#utility for converting text files to .rda to compress output from rjmcmc
#author: JM EASTMAN 2011

cleanup.files <-
function(parmBase, model, fileBase){
	if(model=="BM") {
		parms=c("rates","rate.shifts")
	} else if(model=="OU") {
		parms=c("optima","optima.shifts")
	}
	
	posteriorsamples=lapply(parms, function(x) {
							y=read.table(paste(parmBase,paste(x,"txt",sep="."),sep="/"))
							names(y)=y[1,]
							y=y[-1,]
							row.names(y)=1:nrow(y)
							return(y)
							})
	
	names(posteriorsamples)=parms
	save(posteriorsamples, file=paste(parmBase,paste(fileBase,"posteriorsamples","rda",sep="."),sep="/"))
	unlink(paste(parmBase,paste(parms,"txt",sep="."),sep="/"))
}

