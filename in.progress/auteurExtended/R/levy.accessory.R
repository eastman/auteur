tree.slide=function(phy, node=NULL, up=NA){
	root=Ntip(phy)+1
	edges=cumsum(phy$edge.length)
	ee=edges/max(edges)
	
	# choose random node
	if(is.null(node)){
		choice=phy$edge[min(which(runif(1)<ee)),2]
		choice	
		
	# slide node to one of nearest relative(s) -- based on edge length
	} else if(!is.null(node)){
		if(node==root)stop("Please specify a non-root node.")
		if(is.na(up)) up=as.logical(round(runif(1)))
		if(up){
			choice.tmp=get.ancestor.of.node(node,phy)
			if(choice.tmp==root){
				sisters.tmp=get.desc.of.node(choice.tmp,phy)
				sisters=sisters.tmp[which(sisters.tmp!=node)]
				sister.edges=phy$edge.length[match(sisters, phy$edge[,2])]
				see=cumsum(sister.edges)
				see=see/max(see)
				choice=sisters[min(which(runif(1)<see))]
			} else {
				choice=choice.tmp
			}
		} else {
			descendants=get.desc.of.node(node,phy)
			if(!length(descendants)) {
				choice=tree.slide(phy, node, up=TRUE)
			} else {
				desc.edges=phy$edge.length[match(descendants, phy$edge[,2])]
				dee=cumsum(desc.edges)
				dee=dee/max(dee)
				choice=descendants[min(which(runif(1)<dee))]
			}
		}
	} else {
		stop("If randomly choosing a node, use node=NULL.")
	}	
	return(choice)	
}

jumps.edgewise=function(phy, jump.table){
	jj=rep(0,nrow(phy$edge))
	if(length(jump.table)) {
		jj.tmp=table(jump.table)
		jj[match(names(jj.tmp),phy$edge[,2])]=jj.tmp		
	}
	jj
}



updateJvcv=function(jVCV, from=NULL, to=NULL, node.des){
# from: a node for which a jump has been removed
# to: a node for which a jump has been added
	if(!is.null(from)) {
		ff=node.des[[which(names(node.des)==from)]]
		jVCV[ff,ff]=jVCV[ff,ff]-1
	}
	
	if(!is.null(to)){
		tt=node.des[[which(names(node.des)==to)]]
		jVCV[tt,tt]=jVCV[tt,tt]+1
	}
	jVCV
}

updateLvcv=function(Tvcv, Lvcv, sigsqB, sigsqJ){
	return(sigsqB*Tvcv+sigsqJ*Lvcv)
}

