#utility for generating a list of indicator variables and values for each branch in a phylogeny, using a random model complexity as a starting point (or by constraining model complexity)
#author: JM EASTMAN 2010

generate.starting.point <-
function(data, phy, node.des, model, K=FALSE, jumpsize) { 
	nn=length(phy$edge.length)
	ntip=Ntip(phy)
	if(!K) nshifts=rtpois(1,log(ntip),nn) else nshifts=K-1
	if(nshifts!=0) bb=sample(c(rep(1,nshifts),rep(0,nn-nshifts)),replace=FALSE) else bb=rep(0,nn)
	names(bb)=phy$edge[,2]
	shifts=as.numeric(names(bb)[bb==1])
	if(model=="OU") {
		min.max=c(adjustvalue(min(data), jumpsize), adjustvalue(max(data), jumpsize))
		values=runif(sum(bb)+1, min=min.max[min.max==min(min.max)], max=min.max[min.max==max(min.max)])
	} else if(model=="BM") {
		init.rate=fit.continuous(phy,data)
		min.max=c(adjustrate(init.rate, jumpsize), adjustrate(init.rate, jumpsize))
		values=runif(sum(bb)+1, min=min.max[min.max==min(min.max)], max=min.max[min.max==max(min.max)])
	}
	internal.shifts<-tip.shifts<-numeric(0)
	internal.shifts=sort(shifts[shifts>ntip])
	tip.shifts=shifts[shifts<=ntip]
	
	if(length(internal.shifts)==0 & length(tip.shifts)==0) {
		vv=rep(values, nn)
	} else {
		# assign randomly chosen values for each branch, respecting phylogenetic relationship
		vv=bb
		vv[]=values[length(values)]
		i=0
		if(length(internal.shifts)!=0) {
			for(i in 1:length(internal.shifts)) {
				d=node.des[which(names(node.des)==internal.shifts[i])]
				vv[match(c(internal.shifts[i], unlist(d)), names(vv))]=values[i]
			}
		}
		if(length(tip.shifts)!=0) {
			for(j in 1:length(tip.shifts)) {
				vv[match(tip.shifts[j], names(vv))]=values[j+i]
			}
		}
	}
	return(list(delta=unname(bb), values=unname(vv)))
}

