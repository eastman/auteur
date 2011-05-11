oumat<-function (vmat, alpha, sigma) 
{
	tips=unique(dim(vmat))
	
	out <- .Call("oumat",	
				 TIMES=list(vmat=as.double(array(vmat))),
				 ALPHA=as.double(alpha^2),
				 SIGMA=as.double(sigma), 
				 TIPS=as.integer(tips),
				 PACKAGE = "auteur.extended")
	return(matrix(out$oumat,ncol=tips))
}

		
