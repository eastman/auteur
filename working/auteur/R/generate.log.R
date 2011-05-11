#logging utility used by rjmcmc
#author: JM EASTMAN 2010

generate.log <-
function(bundled.parms, cur.model, file, init=FALSE) {
#	bundled.parms=list(gen=i, mrate=mean(cur.rates), cats=sum(cur.delta.rates)+1, root=cur.root, alpha=cur.alpha, reg=sum(cur.delta.theta)+1, theta=mean(cur.theta), lnL=curr.lnL)
	res=bundled.parms

	if(init) {
		if(cur.model=="BM") msg=paste("gen", "min.rate", "max.rate", "mean.rate", "rates", "root", "lnL", sep="\t") 
	} else {
		if(cur.model=="BM") {
			rr=range(res$rates)
			msg<-paste(res$gen, sprintf("%.3f", rr[1]), sprintf("%.3f", rr[2]), sprintf("%.3f", res$mrate), res$cats, sprintf("%.3f", res$root), sprintf("%.3f", res$lnL),sep="\t")
		} 
	}
	write(msg, file=file, append=TRUE)
}

