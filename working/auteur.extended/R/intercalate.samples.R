#general utility for combining supplied list of dataframes into a single 'intercalated' sample. E.g., list(as.data.frame(c(AAA)),as.data.frame(c(BBB))) becomes data.frame(c(ABABAB))
#author: JM EASTMAN 2010

intercalate.samples <-
function(list.obj) {
	if(length(list.obj)<2) warning("Nothing to be done... too few objects in list.")
	vv=sapply(list.obj, is.vector)
	if(all(vv)) list.obj=lapply(list.obj, function(x) {y=data.frame(x); names(y)=NULL; return(y)})
	orig.names=names(list.obj[[1]])
	n=length(list.obj)
	R=sapply(list.obj, nrow)
	if(length(unique(R))!=1) {
		minR=nrow(list.obj[[which(R==min(R))]])
		list.obj=lapply(list.obj, function(x) {y=as.data.frame(x[1:minR,]); names(y)=orig.names; return(y)})
		warning("Objects pruned to the length of the smallest.")
	}
	C=sapply(list.obj, ncol)
	if(length(unique(C))!=1) stop("Cannot process objects; column lengths do not match.")
	c=unique(C)
	if(c>1) {
		M=lapply(1:length(list.obj), function(x) mm=match(names(list.obj[[x]]), orig.names))
		for(mm in 2:length(list.obj)) {
			list.obj[[mm]]=list.obj[[mm]][,M[[mm]]]
		}
	}
	r=nrow(list.obj[[1]])
	indices=lapply(1:n, function(x) seq(x, r*n, by=n))
	
	out.array=array(dim=c(r*n, c))
	for(nn in 1:n){
		out.array[indices[[nn]],]=as.matrix(list.obj[[nn]])
	}
	
	out.array=data.frame(out.array)
	if(!is.null(orig.names)) names(out.array)=orig.names else names(out.array)=paste("X",1:ncol(out.array),sep="")
	return(out.array)
}

