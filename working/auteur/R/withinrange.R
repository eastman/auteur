#determines if a value falls within a range [min,max]
#author: JM Eastman 2010

withinrange <-
function(x, min, max) {			 
	a=sign(x-min)
	b=sign(x-max)
	if(abs(a+b)==2) return(FALSE) else return(TRUE)
}

