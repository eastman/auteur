#rjmcmc proposal mechanism: scale a rate (limited to only positive values) given a proposal width ('jumpsize')
#author: JM EASTMAN 2010

adjustrate <-
function(rate, jumpsize) {
# mod 10.20.2010 JM Eastman
	vv=rate
	v=log(vv)
	rch <- runif(1, min = -jumpsize, max = jumpsize)
	v=exp(v[]+rch)
	v
}

