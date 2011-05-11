get.ancestors<-function(phy) 
{
	return(lapply(1:Ntip(phy), function(x) get.ancestors.of.node(x,phy)))
}

