#rjmcmc utility for initiating a proposal width for Markov sampling
#author: JM EASTMAN 2010

calibrate.jumpsize <-
function(phy, dat, nsteps=100, model, jumpsizes=NULL) {
	if(!withinrange(nsteps, 100, 1000)) {
		if(nsteps>1000) nsteps=1000 else nsteps=100
	}
	if(is.null(jumpsizes)) jumpsizes=2^(-2:3)
	
	if(model=="BM") {
		acceptance.rates=sapply(jumpsizes, function(x) rjmcmc.bm(phy=phy, dat=dat, ngen=nsteps, jumpsize=x, summary=FALSE, fileBase="jumpsize")$acceptance.rate)
	} else if(model=="OU") {
		acceptance.rates=sapply(jumpsizes, function(x) rjmcmc.ou(phy=phy, dat=dat, ngen=nsteps, jumpsize=x, summary=FALSE, fileBase="jumpsize")$acceptance.rate)
	}
	
	return(sum(jumpsizes*(acceptance.rates/sum(acceptance.rates))))
}

