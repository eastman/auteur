#general statistical function for computing a highest density region (whose proportion is given by 'hpd')
#author: R FITZ-JOHN 2010 (from diversitree)

hdr <-
function(z, hpd, lim=NULL){
	xx <- sort(c(lim, seq(min(z), max(z), length = 1024)))
	ez <- ecdf(z)
	f <- suppressWarnings(approxfun(ez(xx), xx))
	fit <- suppressWarnings(optimize(function(x) f(x + hpd) - f(x), c(0, 1 - hpd)))
	if (inherits(fit, "try-error") || is.na(fit$objective)) stop("HDR interval failure")
	ci <- fit$min
	f(c(ci, ci + hpd))
}


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


#general statistical function for finding probability of a value within a truncated Poisson distribution
#author: JM EASTMAN 2010

ptpois <-
function(x, lambda, k) {
	p.k=ppois(k, lambda, lower.tail=FALSE)
	p.x=ppois(x, lambda, lower.tail=FALSE)
	ptp=p.x/(1-p.k)
	return(ptp)
}


#general statistical function for random samples from a truncated Poisson distribution
#author: JM EASTMAN 2010

rtpois <-
function(N, lambda, k) {
	p=ppois(k, lambda, lower.tail=FALSE)
	out=qpois(runif(N, min=0, max=1-p), lambda)
	out
}



#determines if a value falls within a range [min,max]
#author: JM Eastman 2010

withinrange <-
function(x, min, max) {			 
	a=sign(x-min)
	b=sign(x-max)
	if(abs(a+b)==2) return(FALSE) else return(TRUE)
}



#general statistical function for finding a 'yy' value (by extrapolation or interpolation) given a linear model constructed from 'x' and 'y' and a specified value 'xx' at which to compute 'yy'
#author: JM EASTMAN 2010

infer.y <-
function(xx, x, y) {
	f=lm(y~x)
	p=unname(coef(f))
	yy=xx*p[2]+p[1]
	return(yy)
}

#general utility for finding the indices of the two most similar values in 'x', given 'val' 
#author: JM EASTMAN 2010

nearest.pair <-
function(val, x) {
	dev=abs(x-val)
	names(dev)=1:length(dev)
	return(sort(as.numeric(names(sort(dev))[1:2])))
}

