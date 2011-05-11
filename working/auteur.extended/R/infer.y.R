#general statistical function for finding a 'yy' value (by extrapolation or interpolation) given a linear model constructed from 'x' and 'y' and a specified value 'xx' at which to compute 'yy'
#author: JM EASTMAN 2010

infer.y <-
function(xx, x, y) {
	f=lm(y~x)
	p=unname(coef(f))
	yy=xx*p[2]+p[1]
	return(yy)
}

