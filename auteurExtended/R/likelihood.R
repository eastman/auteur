#primary function for computing the likelihood of data, given a root state, VCV matrix, and Brownian motion model 
#author: LJ HARMON 2009 and JM EASTMAN 2010

bm.lik.fn <-
function(root, orig.dat, vcv, SE) { # from LJ HARMON	
# mod 12.02.2010 JM Eastman: using determinant()$modulus rather than det() to stay in log-space
	y=orig.dat
	
	b <- vcv
	if(any(SE!=0)) diag(b)=diag(b)+SE^2
	w <- rep(root, nrow(b))
	num <- -t(y-w)%*%solve(b)%*%(y-w)/2
	den <- 0.5*(length(y)*log(2*pi) + as.numeric(determinant(b)$modulus))
	list(
		 root=root,
		 lnL = (num-den)[1,1]
		 )	
}


ou.lik.fn <- function(dat, phy, epochs, phyobject, ancestors, vcvSE, rate, alpha, regimes) {
# mod 02.22.2010 JM Eastman
	phyobject$regimes=regimes
	ww <- ouweights(phy, epochs, phyobject, ancestors, alpha)
	w <- ww$W
	v <- oumat(vcvSE, alpha, rate)
	y <- matrix(ww$optima, ncol=1)
	e <- w%*%y - dat
	dim(e) <- n <- length(dat)
	q <- e%*%solve(v,e)
	det.v <- determinant(v,logarithm=TRUE)
	dev <- n*log(2*pi)+as.numeric(det.v$modulus)+q[1,1]
	list(
		 weight=w,
		 vcov=v,
		 resid=e,
		 predict=y,
		 lnL=-0.5*dev
		 )
}


#general phylogenetic utility for determining whether to accept ('r'=1) or reject ('r'=0) a proposed model in Markov chain Monte Carlo sampling
#author: JM EASTMAN 2010

assess.lnR <-
function(lnR) {
	if(is.na(lnR) || abs(lnR)==Inf) {
		error=TRUE
		r=0
	} else {
		if(lnR < -20) {
			r=0 
		} else if(lnR > 0) {
			r=1 
		} else {
			r=exp(lnR)
		}
		error=FALSE
	}
	return(list(r=r, error=error))
}

