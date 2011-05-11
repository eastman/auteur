assignroot <- function(cur.theta, phy) {
	root=phy$edge[1,1]
	desc=phy$edge[which(phy$edge[,1]==root),2]
	root.theta=cur.theta[sample(which(phy$edge[,2]==desc),1)]
	root.theta
}
