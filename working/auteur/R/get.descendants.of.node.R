#general phylogenetic utility for returning the first (usually, unless a polytomy exists) two descendants the supplied node 
#author: JM EASTMAN 2010

get.desc.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}


#general phylogenetic utility for returning all descendants (listed as given in phy$edge[,2]) of a node (excluding the supplied node) 
#author: JM EASTMAN 2010

get.descendants.of.node <-
function(node, phy, tips=FALSE) {
	storage=list()
	if(is.null(node)) node=phy$edge[1,1]
	if(node%in%phy$edge[,1]) {
		desc=c(sapply(node, function(x) {return(phy$edge[which(phy$edge[,1]==x),2])}))
		while(any(desc%in%phy$edge[,1])) {
			desc.true=desc[which(desc%in%phy$edge[,1])]
			desc.desc=lapply(desc.true, function(x) return(get.desc.of.node(x, phy)))
			
			storage=list(storage,union(desc,desc.desc))
			
			desc=unlist(desc.desc)
		}
		storage=as.numeric(sort(unique(unlist(union(desc,storage)))))
	}
	if(tips) return(storage[storage<=Ntip(phy)]) else return(storage)
}

