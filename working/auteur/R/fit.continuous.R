#general phylogenetic utility for rapid computation of the maximum-likelihood estimate of the Brownian motion rate parameter (unless there are polytomies, in which case the function wraps geiger:::fitContinuous)
#author: L REVELL 2010, LJ HARMON 2009, and JM EASTMAN 2010

fit.continuous <-
function(phy, dat){
	n=Ntip(phy)
	ic=try(pic(dat,phy),silent=TRUE)
	if(!inherits(ic, "try-error")) {
		r=mean(ic^2)*((n-1)/n)
		return(r)
	} else {
		return(fitContinuous(phy,dat)$Trait1$beta)
	}
}

