ouch.ou.lik.fn<-function (tree, alpha, sigma, dat, regimes, optima) 
{
#	tree=tree; alpha=alpha; sigma=as.matrix(rate); dat=otd[c("trait")]; regimes=r; optima=unique(r)  
	
	nm <- rownames(dat)
	data <- lapply(as.list(dat), function(x) {
				   names(x) <- nm
				   x
				   })
	
	data <- lapply(data, function(x) x[tree@nodes])
    dat <- do.call(c, lapply(data, function(y) y[tree@term]))	
	
	regimes <- as.factor(regimes)
	names(regimes) <- tree@nodes
	beta <- ouch:::regime.spec(tree, list(regimes))
    n <- length(dat)
    ev <- eigen(alpha, symmetric = TRUE)
    w <- .Call("ouch_weights", object = tree, lambda = ev$values, 
			   S = ev$vectors, beta = beta, PACKAGE="ouch")
    v <- .Call("ouch_covar", object = tree, lambda = ev$values, 
			   S = ev$vectors, sigma.sq = sigma, PACKAGE="ouch")
	y <- matrix(optima, ncol = 1)
    e <- w %*% y - dat
	dim(e) <- n
	
	q <- e %*% solve(v, e)
	det.v <- determinant(v, logarithm = TRUE)
	if (det.v$sign != 1) stop("ou.lik.fn error: non-positive determinant", call. = FALSE)
	dev <- n * log(2 * pi) + as.numeric(det.v$modulus) + q[1, 1]	
    
    list(deviance = dev, coeff = optima, weight = w, vcov = v, 
		 resids = e, lnL = -0.5*dev, beta=beta, ev=ev)
}
