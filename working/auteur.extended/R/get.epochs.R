get.epochs<-function(phy) 
# ordered: tip to root
{
	ancestors=lapply(1:Ntip(phy), function(x) c(x,rev(sort(get.ancestors.of.node(x,phy)))))
	phylo=get.times(phy)
	tt=phylo$time[order(phylo$des)]
	e=lapply(ancestors, function(x) tt[x])
	names(e)=phy$tip.label
	return(e)
}

