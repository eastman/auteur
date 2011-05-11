prepare.ouchdata <- function(phy, data, SE) {
	d=prepare.data(phy, data, SE)
	ape.tre=d$ape.tre
	orig.dat=d$orig.dat
	SE=d$SE
	ot <- ape.to.ouchtree(ape.tre, scale=FALSE)		
	ouch.tre=ot$ouch.tre
	ouch.dat=data.frame(as.numeric(ouch.tre@nodes), as.character(ouch.tre@nodelabels))
	names(ouch.dat)=c("nodes","labels")
	ouch.dat[ouch.dat[,2]=="",2]=NA
	ouch.dat$trait=NA
	mM=match(ouch.dat$labels,names(orig.dat))
	for(i in 1:nrow(ouch.dat)){
		if(!is.na(mM[i])){
			ouch.dat$trait[i]=orig.dat[mM[i]]	
		}
	}
	return(list(ouch.tre=ouch.tre, ouch.dat=ouch.dat, ape.tre=ape.tre, orig.dat=orig.dat, ape.ordered=ot$ape.ordered))
}
