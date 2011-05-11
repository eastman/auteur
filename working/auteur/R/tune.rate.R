#rjmcmc proposal mechanism for updating lineage-specific relative rates (i.e., a subvector swap)
#author: JM EASTMAN 2010

tune.rate <-
function(rates, jumpsize) {
	ss=sample(rates, 1)
	ww=which(rates==ss)
	nn=adjustrate(ss,jumpsize)
	rates[ww]=nn
	return(rates)
}

