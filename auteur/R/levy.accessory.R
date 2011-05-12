tree.slide.LEVY=function(phy, node=NULL, up=NA){
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
		jj=rep(0,nrow(phy$edge))
		jj[match(names(jj.tmp),phy$edge[,2])]=jj.tmp		
	}
	jj
}

adjustjump=function(jump.table, jVCV, k, r, add=FALSE, drop=FALSE, swap=FALSE, phy, node.des){

	if(add==TRUE){
		
		jj=tree.slide(phy)
		jump.table=c(jump.table, jj)
		jVCV=updateJvcv(jVCV, from=NULL, to=jj, node.des)
		todo="add"
	} else if(drop==TRUE | swap==TRUE) {
		
		possible.locs=seq(1,length(jump.table))
		adjustedjump.index<-aji<-possible.locs[ceiling(runif(1)*length(possible.locs))]
		adjustedjump.loc<-ajl<-jump.table[adjustedjump.index]
		
		if(drop==TRUE) {
			
			jump.table=jump.table[-adjustedjump.index]
			jVCV=updateJvcv(jVCV, from=adjustedjump.loc, to=NULL, node.des)
			todo="drop"
		} else if(swap==TRUE){
			
			adjustedjump.newloc<-ajnl<-tree.slide(phy, adjustedjump.loc, up=NA)
			jump.table[adjustedjump.index]=adjustedjump.newloc		
			if((from<-adjustedjump.loc) != (to<-adjustedjump.newloc)) {
				jVCV=updateJvcv(jVCV, from=from, to=to, node.des)
			} else {
				jVCV=jVCV
			}
			todo="swap"
		} else {
			stop("Error encountered.")
		}
	} else {
		stop("Must 'add', 'drop', or 'swap' jumps.")
	}
	
	return(list(jump.table=jump.table, jVCV=jVCV))
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

