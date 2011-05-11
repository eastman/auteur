#general phylogenetic utility for finding the phy$tip.labels associated with the most exclusive monophyletic grouping containing at least supplied 'tips' where 'tips' are elements of the phy$tip.label vector
#author: JM EASTMAN 2010

find.relatives <-
function(tips,phy){
	nn=collect.nodes(tips, phy, strict=FALSE)
	return(phy$tip.label[nn[nn<=Ntip(phy)]])
}

