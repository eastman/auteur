#general phylogenetic plotting utility, given a named vector or data.frame of values that can be associated with phy$edge[,2]
#author: JM EASTMAN 2010
#note: small values are given bluish hues, large values reddish hues; median values are given gray hues

branchcol.plot <-
function(phy, cur.rates, color.length=8, digits=3, plot=TRUE, legend=TRUE, legend.title="", ...) {
	if(any(cur.rates<=0)) log=FALSE else log=TRUE
	require(ape)
	if(!is.null(names(cur.rates))) {
		cur.rates=cur.rates[match(as.numeric(names(cur.rates)), phy$edge[,2])] 
	} else {
		names(cur.rates)=phy$edge[,2]
		warning("Rates assumed to be ordered as in 'phy$edge'")
	}
	cur.rates=sapply(cur.rates, mean)
	require(colorspace)
	if(log) ests=log(cur.rates) else ests=cur.rates
	mm=sapply(ests, function(x) x-median(ests))
	cce=diverge_hcl(2*color.length+1, power = 0.5)
	e.seq=seq(-max(abs(mm)),max(abs(mm)),length=2*color.length+1)
	colors.branches=sapply(mm, function(x) cce[which(min(abs(e.seq-x))==abs(e.seq-x))])
	lseq=e.seq+median(ests)
	lseq=seq(min(lseq), max(lseq), length=color.length)
	lcce=cce[round(seq(1, length(cce), length=color.length))]
	if(log) lseq=exp(rev(lseq)) else lseq=rev(lseq)
	
	if(plot) {
		plot.phylo(phy, cex=0.1, edge.color=colors.branches, ...)
		if(legend) {
			legend("topright", title=legend.title, cex=0.5, pt.cex=1, text.col="darkgray", 
				   legend = sprintf(paste("%", 2*digits, paste(digits, "f", sep=""), sep="."), lseq), pch=21, ncol=1, col = "darkgray", 
				   pt.bg = rev(lcce), box.lty="blank", border="white")
		}
	} else {
		return(list(col=colors.branches,legend.seq=lseq,legend.col=lcce))
	}
}

