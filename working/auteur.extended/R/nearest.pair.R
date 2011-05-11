#general utility for finding the indices of the two most similar values in 'x', given 'val' 
#author: JM EASTMAN 2010

nearest.pair <-
function(val, x) {
	dev=abs(x-val)
	names(dev)=1:length(dev)
	return(sort(as.numeric(names(sort(dev))[1:2])))
}

