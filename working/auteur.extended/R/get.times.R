get.times<-function(phy) 
{
	df=data.frame(phy$edge)
	names(df)=c("anc","des")
	anc=lapply(df$des,get.ancestors.of.node,phy)
	order.edges=c(0,phy$edge.length)[order(c(Ntip(phy)+1,df$des))]
	df$time=0
	for(n in 1:length(anc)){
		df$time[n]=sum(order.edges[anc[[n]]])+order.edges[df$des[n]]
	}
	
	df=rbind(c(0,Ntip(phy)+1,0),df)

	return(df)
}

