#utility within rjmcmc for determining which proposal was accepted in a Markov generation, given a vector of proposal acceptances
#author: JM EASTMAN 2010

determine.accepted.proposal <-
function(startparms, endparms, cOK) {
	which(startparms!=endparms)->marker
	cOK[marker]=cOK[marker]+1
	cOK
}

