		
ouweights <- function(phy, epochs, object, ancestors, alpha){
	if(!all(c("anc","des","time","regimes")%in%names(object)))stop("Cannot decipher supplied 'object'.")
	ntip = Ntip(phy)
	object=object[order(object$des),]
	optima = sort(unique(as.numeric(object$regimes)))
    nreg = length(optima)
	
	out <- .Call("ouweights",	
				 PHYLO=list(
							epochs=epochs,
							ancestors=ancestors,
							regimes=as.double(object$regimes),
							optima=as.double(optima),
							alpha=as.double(alpha^2),
							tips=as.integer(ntip)),
				 PACKAGE = "auteur.extended")
	
	W=matrix(out$weights,ncol=nreg)
	opt=optima
	return(list(W=W, optima=opt))
} 