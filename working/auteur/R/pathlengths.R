pathlengths=function(phy){
tips=1:Ntip(phy)
root=Ntip(phy)+1
addlength=function(length, phy, node){
	a=get.ancestor.of.node(node, phy)
	if(a!=root) {
		addlength(length+phy$edge.length[which(phy$edge[,2]==a)], phy, a)
	} else {
		return(length)
	}
	}
paths=sapply(tips, function(x) addlength(phy$edge.length[which(phy$edge[,2]==x)],phy,x))
return(paths)	
}
