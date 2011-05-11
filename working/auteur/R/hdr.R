#general statistical function for computing a highest density region (whose proportion is given by 'hpd')
#author: R FITZ-JOHN 2010 (from diversitree)

hdr <-
function(z, hpd, lim=NULL){
	xx <- sort(c(lim, seq(min(z), max(z), length = 1024)))
	ez <- ecdf(z)
	f <- suppressWarnings(approxfun(ez(xx), xx))
	fit <- suppressWarnings(optimize(function(x) f(x + hpd) - f(x), c(0, 1 - hpd)))
	if (inherits(fit, "try-error") || is.na(fit$objective)) stop("HDR interval failure")
	ci <- fit$min
	f(c(ci, ci + hpd))
}

