#intercalates samples from multiple independent Markov chains (generated by rjmcmc.bm())
#author: JM EASTMAN 2010

pool.rjmcmcsamples <-
function(base.dirs, lab="combined"){
	outdir=paste(lab,"combined.rjmcmc",sep=".")
	if(!file.exists(outdir)) dir.create(outdir)
	
	# estimated parameters
	samples=lapply(base.dirs, function(x) {load(paste(x, dir(x, pattern="posteriorsamples.rda"), sep="/")); return(posteriorsamples)})
	nn=unique(unlist(lapply(samples, names)))
	posteriorsamples=list()
	for(i in nn) {
		res.sub=intercalate.samples(lapply(1:length(samples), function(x) samples[[x]][[i]]))
		posteriorsamples[[which(nn==i)]]=res.sub
	}
	names(posteriorsamples)=nn
	save(posteriorsamples, file=paste(outdir, paste(lab,"posteriorsamples.rda",sep="."), sep="/"))
	
	# logfile
	
	write.table(intercalate.samples(lapply(base.dirs, function(x) {y=read.table(paste(x, dir(x, pattern="rjmcmc.log"), sep="/"), header=TRUE); return(y)})), paste(outdir, paste(lab,"rjmcmc.log",sep="."), sep="/"), quote=FALSE, row.names=FALSE)
}

