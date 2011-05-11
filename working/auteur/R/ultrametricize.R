
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