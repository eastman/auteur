#general statistical function for random samples from a truncated Poisson distribution
#author: JM EASTMAN 2010

rtpois <-
function(N, lambda, k) {
	p=ppois(k, lambda, lower.tail=FALSE)
	out=qpois(runif(N, min=0, max=1-p), lambda)
	out
}

