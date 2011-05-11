#general phylogentic utility wrapping geiger:::treedata for checking data consistency
#author: JM EASTMAN 2010

prepare.data <-
function(phy, data, SE) {
	td <- treedata(phy, data, sort = TRUE)	
	if(any(SE!=0)) {
		if(!is.null(names(SE))) {
			se=rep(0,length(td$phy$tip.label))
			names(se)=td$phy$tip.label
			ss.match=match(names(SE), names(se))
			if(any(is.na(ss.match))) warning(paste(names(SE[is.na(ss.match)]), "not found in the dataset", sep=" "))
			se[match(names(SE[!is.na(ss.match)]),names(se))]=SE[!is.na(ss.match)]
			SE=se
		} else {
			if(name.check(phy,data)=="OK" & length(data)==length(SE)) {
				names(SE)=td$phy$tip.label
				warning("ordering of values in SE assumed to be in perfect correspondence with phy$tip.label")
			} else {
				stop("SE must be a named vector of measurement error")
			}
		}
	} else {
		SE=rep(0,length(td$phy$tip.label))
		names(SE)=td$phy$tip.label
	}
	
	return(list(ape.tre=td$phy, orig.dat=td$data[,1], SE=SE))
}

