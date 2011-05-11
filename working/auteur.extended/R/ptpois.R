#general statistical function for finding probability of a value within a truncated Poisson distribution
#author: JM EASTMAN 2010

ptpois <-
function(x, lambda, k) {
	p.k=ppois(k, lambda, lower.tail=FALSE)
	p.x=ppois(x, lambda, lower.tail=FALSE)
	ptp=p.x/(1-p.k)
	return(ptp)
}

