#rjmcmc proposal mechanism: scale a value given a proposal width ('jumpsize')
#author: JM EASTMAN 2010

adjustvalue <-
function(value, jumpsize, theta=FALSE) {
# mod 10.20.2010 JM Eastman
	vv=value
	if(runif(1)<0.5 | length(unique(vv))==1) {
		rch <- runif(1, min = -jumpsize, max = jumpsize)
		return(vv[]+rch)
	} else {
		return((vv-mean(vv))*runif(1,min=-jumpsize, max=jumpsize)+mean(vv))
	}
}

