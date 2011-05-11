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
