#general statistical function for comparing two vectors on non-independent samples
#author: JM EASTMAN 2010

randomization.test <-
function(obs=obs, exp=exp,  mu=0, iter=10000, two.tailed=FALSE){
	
	O=as.numeric(obs)[sample(1:length(obs), iter, replace=TRUE)]
	E=as.numeric(exp)[sample(1:length(exp), iter, replace=TRUE)]
	
	result=c(O-(E+mu))
	p=round(sum(result>=0)/iter, digits=nchar(iter))
	q=1-p
	
	if(two.tailed) {
		res=list(diffs=result, 2*(c(p,q)[which(c(p,q)==min(c(p,q)))]))
	} else {
		res=list(diffs=result, greater=p, lesser=q)
	}
	return(res)
}

