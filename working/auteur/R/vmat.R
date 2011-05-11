vmat <- function(phy){
	n=Ntip(phy)
	out <- .Call("vmat", tree=list(
									ROOT = as.integer(n+1),
									MAXNODE = as.integer(max(phy$edge[,1])),
									ENDOFCLADE = as.integer(dim(phy$edge)[1]),
									ANC = as.integer(phy$edge[,1]),
									DES = as.integer(phy$edge[,2]),
									EDGES = as.double(c(phy$edge.length,0)),
									VCV = as.double(array(matrix(0, n, n)))),
									PACKAGE = "auteur")
	v=matrix(out$VCV,nrow=n,byrow=FALSE)
	rownames(v)<-colnames(v)<-phy$tip.label
	return(v)
} 

		
