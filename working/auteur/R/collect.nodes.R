#general phylogenetic utility for finding the edge labels (from phy$edge[,2]) associated with the most exclusive monophyletic grouping containing at least supplied 'tips' where 'tips' are elements of the phy$tip.label vector
#author: JM EASTMAN 2010

collect.nodes <-
function(tips,phy,strict=FALSE){
	if(length(tips)==1) stop("Must supply at least two tips.")
	tt=phy$tip.label
	tt=lapply(tips, function(x) which(phy$tip.label==x))
	nodes=lapply(tt, function(x) get.ancestors.of.node(x, phy))
	all.nodes=nodes[[1]]
	for(i in 2:length(nodes)) {
		all.nodes=intersect(all.nodes,nodes[[i]])
	}
	base.node=max(all.nodes)
	if(strict) {
		u=unlist(nodes)
		return(c(sort(unique(u[u>=base.node])), unlist(tt)))
	} else {
		return(c(base.node,get.descendants.of.node(base.node,phy)))
	}
}

