#general phylogenetic utility for determining whether a node is the root of the phylogeny
#author: JM EASTMAN 2010

is.root <-
function(node,phy) {
	if(node==phy$edge[1,1]) return(TRUE) else return(FALSE)
}

