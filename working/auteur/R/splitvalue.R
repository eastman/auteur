#rjmcmc proposal mechanism: split a value into two while maintaining the mean
#author: JM EASTMAN 2010

splitvalue <-
function(cur.vv, n.desc, n.split, factor=log(2)){
	dev=cur.vv-adjustvalue(cur.vv, factor)
	nr.desc=cur.vv + dev/n.desc
	nr.split=cur.vv - dev/n.split
	return(list(nr.desc=nr.desc, nr.split=nr.split))
}

