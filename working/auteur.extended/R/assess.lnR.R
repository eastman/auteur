#general phylogenetic utility for determining whether to accept ('r'=1) or reject ('r'=0) a proposed model in Markov chain Monte Carlo sampling
#author: JM EASTMAN 2010

assess.lnR <-
function(lnR) {
	if(is.na(lnR) || abs(lnR)==Inf) {
		error=TRUE
		r=0
	} else {
		if(lnR < -20) {
			r=0 
		} else if(lnR > 0) {
			r=1 
		} else {
			r=exp(lnR)
		}
		error=FALSE
	}
	return(list(r=r, error=error))
}

