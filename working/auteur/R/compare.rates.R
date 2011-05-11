#general phylogenetic statistical utility for comparing values ('posterior.rates') from posterior densities for two sets of supplied branches ('nodes')
#author: JM EASTMAN 2010

compare.rates <-
function(branches=list(A=c(NULL,NULL),B=c(NULL,NULL)), phy, posterior.values, ...){
	nodes=branches
	rates=posterior.values
	if(is.null(names(nodes))) names(nodes)=c("A","B")
	rates.a=rates[,match(nodes[[1]],names(rates))]
	rates.b=rates[,match(nodes[[2]],names(rates))]
	edges.a=phy$edge.length[match(nodes[[1]],phy$edge[,2])]
	edges.b=phy$edge.length[match(nodes[[2]],phy$edge[,2])]
	weights.a=edges.a/sum(edges.a)
	weights.b=edges.b/sum(edges.b)
	if(length(nodes[[1]])>1) {r.scl.a=apply(rates.a, 1,function(x) sum(x*weights.a))} else {r.scl.a=sapply(rates.a, function(x) sum(x*weights.a))}
	if(length(nodes[[2]])>1) r.scl.b=apply(rates.b, 1,function(x) sum(x*weights.b)) else r.scl.b=sapply(rates.b, function(x) sum(x*weights.b))
	return(list(scl.A=r.scl.a, scl.B=r.scl.b, r.test=randomization.test(r.scl.a, r.scl.b, ...)))	
}

