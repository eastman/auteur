#general phylogenetic utility for returning first ancestor (as a numeric referent) of the supplied node 
#author: JM EASTMAN 2010

get.ancestor.of.node <-
function(node, phy) {
	anc=phy$edge[which(phy$edge[,2]==node),1]
	anc
}



#general phylogenetic utility for returning all ancestors (listed as given in phy$edge[,2]) of a node 
#author: JM EASTMAN 2010

get.ancestors.of.node <-
function(node, phy) {
	ancestors=c()
	a=0
	while(1) {
		a=a+1
		node<-as.numeric(get.ancestor.of.node(node, phy))
		if(!length(node)) break() else ancestors[a]=node
	}
	return(ancestors)
}

