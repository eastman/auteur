#scales VCV matrix of phylogeny by relative evolutionary rates under Brownian motion
#author: JM EASTMAN 2010

updatevcv <-
function(ape.tre, new.rates) {
	n=ape.tre
	n$edge.length=n$edge.length*new.rates
	vv=vmat(n)
	return(vv)
}

