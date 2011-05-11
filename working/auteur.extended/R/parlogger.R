#logging utility used by rjmcmc
#author: JM EASTMAN 2010

parlogger <-
function(model, init=FALSE, node.des, i, curr.lnL, cur.root=NULL, cur.rates=NULL, cur.delta.rates=NULL, cur.alpha=NULL, cur.theta=NULL, cur.delta.theta=NULL, parmBase) {	
	if(init) {
		msg.long=names(node.des[-1])
		if(model=="BM") {
			parlogs=paste(parmBase, paste(c("rates", "rate.shifts"), "txt", sep="."), sep="")
			sapply(parlogs[1:length(parlogs)], function(x) write.table(paste(msg.long,collapse=" "), x, quote=FALSE, col.names=FALSE, row.names=FALSE))
		} else {
			parlogs=paste(parmBase, paste(c("optima", "optima.shifts"), "txt", sep="."), sep="")
			sapply(parlogs[1:length(parlogs)], function(x) write.table(paste(msg.long,collapse=" "), x, quote=FALSE, col.names=FALSE, row.names=FALSE))
		}
	} else {
		if(model=="BM") {
			parlogs=paste(parmBase, paste(c("rates", "rate.shifts"), "txt", sep="."), sep="")
			outlist=list(cur.rates, cur.delta.rates) 
			sapply(1:length(outlist), function(x) write.table(paste(outlist[[x]],collapse=" "), parlogs[[x]], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)) 
		} else {
			parlogs=paste(parmBase, paste(c("optima", "optima.shifts"), "txt", sep="."), sep="")
			outlist=list(cur.theta, cur.delta.theta) 
			sapply(1:length(outlist), function(x) write.table(paste(outlist[[x]],collapse=" "), parlogs[[x]], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)) 
		}
	}
}

