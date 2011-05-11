#primary function for computing the likelihood of data, given a root state, VCV matrix, and Brownian motion model 
#author: LJ HARMON 2009 and JM EASTMAN 2010

bm.lik.fn <-
function(rates, root, orig.dat, vcv, SE) { # from LJ HARMON	
# mod 12.02.2010 JM Eastman: using determinant()$modulus rather than det() to stay in log-space
	y=orig.dat
	
	b <- vcv
	if(any(SE!=0)) diag(b)=diag(b)+SE^2
	w <- rep(root, nrow(b))
	num <- -t(y-w)%*%solve(b)%*%(y-w)/2
	den <- 0.5*(length(y)*log(2*pi) + as.numeric(determinant(b)$modulus))
	list(
		 lnL = (num-den)[1,1]
		 )	
}

