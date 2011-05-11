#rjmcmc proposal mechanism for updating lineage-specific relative values (i.e., a subvector swap)
#author: JM EASTMAN 2010

tune.value <-
function(values, jumpsize) {
	ss=sample(values, 1)
	ww=which(values==ss)
	nn=adjustvalue(ss,jumpsize)
	values[ww]=nn
	return(values)
}

