#plotting function for comparing posterior densities of estimates
#author: JM EASTMAN 2010

trace.plot <-
function(obj, col, alpha, lwd=1, hpd=0.95, bars=TRUE, legend.control=list(plot=TRUE, pos=NA, cex=1, pt.cex=1, pch=22, title=""), truncate=list(min=NA, max=NA), xlim=list(min=NA, max=NA), ...){
	
	if(!is.data.frame(obj)) {
		if(is.vector(obj)) obj=as.data.frame(obj) else stop("Object must be supplied as a vector or dataframe.")
	}
	
	# prepare (or assign) necessary arguments 
	control=list(plot=TRUE,pos=NA,cex=1,pt.cex=1,pch=22,title="")
	control[names(legend.control)]<-legend.control
	legend.control=control
	if(missing(col)) col=gray.colors(ncol(obj))
	if(length(col)!=ncol(obj)) col=gray.colors(ncol(obj))
	if(missing(alpha)) alpha=0.8
	line.col=add.transparency(col, min(c(0.95, alpha+alpha*0.25)))
	col=add.transparency(col, alpha)
	
	# truncate data
	if(any(!is.na(truncate))) {
		trunc.todo=c("min","max")[match(names(!is.na(truncate)),c("min","max"))]
		for(t in 1:length(trunc.todo)) {
			if(length(truncate[[t]])!=ncol(obj) & length(truncate[[t]]==1)) {
				truncate[[t]]=rep(truncate[[t]],ncol(obj))
				warning(paste("assuming truncation value for ", sQuote(trunc.todo[t]), " is the same for all columns in 'obj'.", sep=""))
			} 
			for(c in 1:ncol(obj)) {
				if(trunc.todo[t]=="min") {
					if(!is.na(truncate$min[c])) {
						mi=sapply(obj[,c], function(x) if(is.na(x) | x<truncate$min[c]) return(TRUE) else return(FALSE))
						obj[mi,c]=NA
					}
				} else if(trunc.todo[t]=="max") {
					if(!is.na(truncate$max[c])) {
						ma=sapply(obj[,c], function(x) if(is.na(x) | x>truncate$max[c]) return(TRUE) else return(FALSE))
						obj[ma,c]=NA
					}
				}
			}
		}
	}
	
	# prepare density objects and compute HDRs 
	dd=lapply(1:ncol(obj), function(x) {y=range(obj[,x],na.rm=TRUE); density(obj[,x], from=min(y), to=max(y), na.rm=TRUE, n=1024, kernel="cosine")})
	if(!is.null(hpd)) {
		hh=lapply(1:ncol(obj), function(z) {zz=obj[!is.na(obj[,z]),z]; hdr(zz, hpd=hpd)})
	} else {
		hh=lapply(1:ncol(obj), function(z) {zz=obj[!is.na(obj[,z]),z]; c(min(zz), max(zz))})
	}
	xx=lapply(1:ncol(obj),function(z) {zz=obj[!is.na(obj[,z]),z]; return(c(min(zz[which(zz>=hh[[z]][1])]), max(zz[which(zz<=hh[[z]][2])])))})
	density.yy=lapply(1:length(dd), function(z) {sapply(xx[[z]], function(w) {xs=nearest.pair(w, dd[[z]]$x); infer.y(w, dd[[z]]$x[xs], dd[[z]]$y[xs])})})
	
	# include HDR points in density objects (necessary for precision in plotting)
	dd.new=lapply(1:length(dd), function(z) {
				  xxx=c(xx[[z]],dd[[z]]$x)
				  yyy=c(density.yy[[z]],dd[[z]]$y)
				  dx=sort(c(xx[[z]],dd[[z]]$x),index=TRUE)
				  return(dd.out=list(x=dx$x, y=yyy[dx$ix]))
				  })
	for(z in 1:length(dd)) {
		dd[[z]]$x=dd.new[[z]]$x
		dd[[z]]$y=dd.new[[z]]$y
	}
	
	ii=lapply(1:length(dd), function(z) sapply(dd[[z]]$x, function(w) withinrange(w, min(hh[[z]]), max(hh[[z]]))))
	
	# find xlims for plotting 
	lims.x.tmp=unlist(lapply(1:length(dd), function(z) return(c(max(dd[[z]]$x), min(dd[[z]]$x)))))
	lims.x=c(min(lims.x.tmp), max(lims.x.tmp))
	lims.x=c(min(lims.x)-0.1*min(lims.x), max(lims.x)+0.1*max(lims.x))
	if(any(!is.na(xlim))) {
		xlim.todo=c("min","max")[match(names(!is.na(xlim)),c("min","max"))]
		for(t in 1:length(xlim.todo)) {
			if(xlim.todo[t]=="min") lims.x[1]=xlim$min
			if(xlim.todo[t]=="max") lims.x[2]=xlim$max
		}
	}
	
	lims.y.tmp=unlist(lapply(1:length(dd), function(z) return(c(max(dd[[z]]$y), min(dd[[z]]$y)))))
	lims.y=c(-0.05, 1.05)*max(lims.y.tmp)
	
	# PLOTTING of densities and legend  
	plot(x=NULL, xlim=lims.x, ylim=lims.y, ylab="density", bty="n", type="n", ...)
	q=seq(0,min(lims.y)+0.4*min(lims.y),length=ncol(obj)+2)
	q=q[-c(1:2)]
	for(i in 1:length(dd)) {
		dat=dd[[i]]
		index=ii[[i]]
		xs=xx[[i]]
		ys=density.yy[[i]]
		polygon(c(dat$x[index], rev(dat$x[index])), c(dat$y[index],rep(0, length(dat$y[index]))), col=col[i], border=NA)
		dat$x=c(min(dat$x, na.rm=TRUE),dat$x,max(dat$x, na.rm=TRUE))
		dat$y=c(0,dat$y,0)
		lines(dat, col=line.col[i], lwd=lwd)
		if(bars) {
			arrows(xx[[i]][1], q[i], xx[[i]][2], q[i], code = 1, length = 0.0, col = line.col[i], lwd=lwd)
			points(data.frame(list(c(xx[[i]][1], xx[[i]][2]), c(q[i], q[i]))), pch=21, col=line.col[i], bg=line.col[i], cex=0.5*lwd)
		}
	}
	if(all(!is.null(names(obj))) & legend.control$plot) {
		if(is.na(legend.control$pos)) {
			mm=c(unlist(obj))
			mm=mm[!is.na(mm)]
			if(abs(mean(mm)-min(mm))>abs(mean(mm)-max(mm))) pos="topleft" else pos="topright"
		} else {
			pos=legend.control$pos
		}
		legend(x=pos, names(obj), pch=legend.control$pch, col="gray", pt.bg=line.col, bty="n", cex=legend.control$cex, title=legend.control$title, pt.cex=legend.control$pt.cex)	
	}
}

