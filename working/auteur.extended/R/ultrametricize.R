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
get.ancestor.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,2]==node),1])
}
ultrametricize=function(phy, tol=1e-8, trim=c("min","max","mean")){
	paths=pathlengths(phy)
	trim=switch(match.arg(trim),
         min = min(paths),
         max = max(paths),
         mean = mean(paths))
	if(diff(range(paths))<=tol) {
		for(i in 1:length(paths)){
			e=phy$edge.length[which(phy$edge[,2]==i)->bl]
			phy$edge.length[bl]=e+(trim-paths[i])
			}
		return(phy)
	} else {
		stop("Difference in path lengths is larger than supplied tolerance")
	}
}