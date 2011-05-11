#rjmcmc run diagnosis, generating plots to compare posterior and prior densities of K-rate models and conduct Bayes factor comparisons between models differing in complexity
#author: JM EASTMAN 2010

tracer <-
function(base.log, lambdaK=log(2), Bayes.factor=NULL, burnin=0.25, col.line=1, pdf=TRUE, factor="M",...) {
	d=read.table(base.log,header=TRUE)
	tag.tmp=strsplit(base.log,".",fixed=TRUE)[[1]]
	tag=paste(tag.tmp[(tag.tmp != "rjTrait" & tag.tmp != "log")],collapse=".")
	
	if(pdf) pdf(paste(tag,"pdf",sep="."))
	
	param=names(d)[5]
	sparam=ifelse(param=="rates","rate","optimum")
	
# lnL trace
	if(any(names(d)=="lnL")){
		d$col=1
		if(length(unique(d$state))!=nrow(d)) {
			cols=nrow(d)/length(unique(d$state))
			if(cols%%1==0) {
				d$col=rep(1:cols, nrow(d)/cols)
			}
		} 
		grays=gray.colors(max(d$col),start=0,end=0.6)
		d$col=grays[d$col]
		if(!is.null(factor)) {
			if(factor=="k" | factor=="K") {
				plot(d$state/1000, d$lnL, main="likelihood trace", ylab="log-likelihood", xlab="generation (k)", bty="n", cex=0.2, col=add.transparency(d$col,0.75), bg=add.transparency(d$col,0.75), pch=21)
			} else {
				plot(d$state/1000000, d$lnL, main="likelihood trace", ylab="log-likelihood", xlab="generation (M)", bty="n", cex=0.2, col=add.transparency(d$col,0.75), bg=add.transparency(d$col,0.75), pch=21)
			}
		} else {
			plot(d$state, d$lnL, main="likelihood trace", ylab="log-likelihood", xlab="generation", bty="n", cex=0.2, col=add.transparency(d$col,0.75), bg=add.transparency(d$col,0.75), pch=21)	
		}
	}
	
# mean density
	if(any(names(d)=="mean")) dmr=as.data.frame(d$mean) 
	if(length(col.line)!=ncol(dmr)) {
		col.line=gray.colors(ncol(dmr)+1)
		col.line=col.line[-length(col.line)]
	}
	if(!pdf) dev.new()
	if(all(dmr>0)) {
		trace.plot(dmr, col=col.line, xlab=param, hpd=NULL, legend.control=list(plot=FALSE), main=paste("posterior sample of mean ", param, sep=""), ...)	
	} else {
		trace.plot(dmr, col=col.line, xlab=param, hpd=NULL, legend.control=list(plot=FALSE), main=paste("posterior sample of mean ", param, sep=""), ...)	
	}
	
# K-parameter model prior and posterior probabilities; Bayes factor comparison
	if(any(names(d)==param)) {
		bb=max(1, round(burnin*nrow(d)))
		d=d[bb:nrow(d),]
		prior=rpois(nrow(d), lambdaK)+1
		tt.prior=table(prior)
		tt.post=table(d[,param])
		df=array(dim=c(2,max(c(prior,d[,param]))))
		for(pr in 1:length(tt.prior)) {
			df[1,(as.numeric(names(tt.prior)[pr]))]=tt.prior[pr]
		}
		for(pp in 1:length(tt.post)) {
			df[2,(as.numeric(names(tt.post)[pp]))]=tt.post[pp]
		}
		df[is.na(df)]=0
		dimnames(df)[[2]]=seq(1,ncol(df))
		if(!pdf) dev.new()
		df=apply(df,2,function(x)x/sum(df/2))
		
# posterior-prior plot of K-rate models
		barplot(df,col=c("whitesmoke","black"),xlab=param,ylab="probability",ylim=c(0,1),beside=TRUE,main=paste("probabilities of k-",param, " models",sep=""))
		legend(x="topright", c("prior","posterior"), fill=c("whitesmoke","black"), bty="n", cex=0.75)
		
# Bayes factor comparison
		if(!is.null(Bayes.factor)) {
			is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
			resample <- function(a=c(0,0), b=vector(), resamples=1000, frequency=0.5, multi=FALSE) {
				if(!multi) {
					xx=lapply(1:resamples, function(x) {
							  rr=sample(b, size=round(frequency*length(b)), replace=TRUE)
							  return(sum(rr==a[2])/sum(rr==a[1]))
							  }
							  )
				} else {
					xx=lapply(1:resamples, function(x) {
							  rr=sample(b, size=round(frequency*length(b)), replace=TRUE)
							  return(sum(rr>1)/sum(rr==1))
							  }
							  )
				}
				return(unlist(xx))
			}
			if(length(Bayes.factor)==2 ) {
				if(!pdf) dev.new()
				if(Bayes.factor[2]!="multi") {
					BF.tmp=sapply(Bayes.factor, function(x)if(sum(d[,param]==x)>0)return(TRUE)else return(FALSE))
					if(all(BF.tmp)) {
						BF=sum(d[,param]==Bayes.factor[2])/sum(d[,param]==Bayes.factor[1])
						resampled.BF=resample(Bayes.factor, d[,param])
						trace.plot(resampled.BF[is.finite(resampled.BF)], legend.control=list(plot=FALSE), main=paste("Bayes factor: ", paste(paste(rev(Bayes.factor), "-", sep=""), collapse=" versus "), paste(sparam,"model of evolution",sep=" "), sep=""), xlab="Bayes factor")
						lines(data.frame(c(BF,BF),c(0,max(density(resampled.BF)$y)), lty=4))
					} else {
						BF=paste("Model(s)", paste(Bayes.factor[which(!BF.tmp)], collapse=" and "), "appear(s) to be missing from the posterior sample.", sep=" ")
						resampled.BF=NA
					}
				} else {
					BF.tmp=c((sum(d[,param]==1)>0), (sum(d[,param]>1)>0))
					if(all(BF.tmp)) {
						BF=sum(d[,param]>1)/sum(d[,param]==1)
						resampled.BF=resample(Bayes.factor, d[,param], multi=TRUE)
						trace.plot(resampled.BF[is.finite(resampled.BF)], legend.control=list(plot=FALSE), main=paste("Bayes factor: ", paste(paste(c("multi","global"), "-", sep=""), collapse=" versus "), paste(sparam,"model of evolution",sep=" "), sep=""), xlab="Bayes factor")
						lines(data.frame(c(BF,BF),c(0,max(density(resampled.BF)$y)), lty=4))
					} else {
						BF=paste("Model(s)", paste(Bayes.factor[which(!BF.tmp)], collapse=" and "), "appear(s) to be missing from the posterior sample.", sep=" ")
						resampled.BF=NA
					}
				}
				
			}
		}
	}
	
	if(pdf) dev.off()	
	if(!is.null(Bayes.factor)) return(list(BF=BF, resampled.BF=resampled.BF))
}



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



#plotting function for comparing posterior densities of evolutionary process 'shifts' along a phylogeny
#author: JM EASTMAN 2010

process.shifts<-function(phy, shifts, level) {
	if(!all(names(shifts)%in%phy$edge[,2])) stop("Cannot process shifts: check that tree matches posterior summaries")
	hits.tmp<-apply(shifts,1,sum)
	hits=length(hits.tmp[hits.tmp>0])
	branches.tmp=apply(shifts, 2, function(x) sum(x))
	branches.tmp=branches.tmp/hits
	branches=branches.tmp[branches.tmp>=level]
	branches=branches[order(branches, decreasing=TRUE)]
	nn=as.numeric(names(branches))
	desc=lapply(nn, function(x) {y=get.descendants.of.node(x, phy=phy, tips=TRUE); if(length(y)) return(phy$tip.label[y]) else return(NULL)})
	names(desc)=nn
	if(length(nn) & any(!is.na(nn))) {
		di=sapply(nn, function(x) if(x>Ntip(phy)) return(TRUE) else return(FALSE))
		if(any(di)) desc=desc[di] else desc=NULL
	} else {
		desc<-branches<-NULL
	}
	return(list(descendants=desc, branch.shift.probability=branches))
}

shifts.plot <-
function(phy, base.dir, burnin=0, level=0.01, internal.only=FALSE, paint.branches=TRUE, legend=TRUE, pdf=TRUE, verbose=TRUE, lab="OUT", ...) {
	color.length=17
	oldwd=getwd()
	setwd(base.dir)
	
	if(pdf | verbose) {
		out.base=paste(lab, paste("burnin", burnin, sep=""), paste("fq", level, sep=""), sep=".")
		outdir=paste("./results")
		if(!file.exists(outdir)) dir.create(outdir)
	}
	
	# collect data
	cat("READING estimates...\n")
	posteriorsamples=NULL
	if(length(dir(,pattern="posteriorsamples.rda"))==1) load(dir()[grep("posteriorsamples.rda",dir())]) else stop("Estimates cannot be located: .rda file appears missing")
	shifts=posteriorsamples$shifts
	burnin=ceiling(burnin*nrow(shifts))
	shifts=shifts[-c(1:(burnin)),]
	shifts.res=process.shifts(phy, shifts, level)
	ests=posteriorsamples[[(which(names(posteriorsamples)!="shifts")->target)]]
	ests=ests[-c(1:(burnin)),]
	average.ests.tmp=apply(ests,2,mean)
	average.ests=average.ests.tmp-median(average.ests.tmp)
	
	# determine whether to use logspace for plotting
	param=names(posteriorsamples)[target]
	if(param%in%c("optima","bursts")) logspace=FALSE else logspace=TRUE
	
	# collect edge colors (for rates)
	if(paint.branches) {
		colors.branches.tmp=branchcol.plot(phy, ests, plot=FALSE, color.length=color.length, log=logspace)
		colors.branches=colors.branches.tmp$col
	} else {
		colors.branches=1
	}
	
	# collect node colors (for shifts)
	ccx=diverge_hcl(color.length, power = 0.5)
	c.seq=round(seq(-1,1,length=color.length),digits=1)
	all.nodes=seq(1:(Ntip(phy)+Nnode(phy)))[-(Ntip(phy)+1)]

	if(!is.null(shifts.res$branch.shift.probability)) {
		require(colorspace)
		nodes=as.numeric(names(shifts.res$branch.shift.probability))
		shift.direction=sapply(all.nodes, function(x) {
					a=get.ancestor.of.node(x, phy)
					comp=ests[,which(names(ests)==a)]
					this=ests[,which(names(ests)==x)]
					if(!length(comp)) { # dealing with first descendant of root
						d=get.desc.of.node(a, phy)
						d=d[which(d!=x)]
						d.shifts=shifts[, which(names(shifts)==d)]
						comp=ests[,which(names(ests)==d)]
						x.res=sapply(1:length(d.shifts), function(y) {
								 if(d.shifts[y]==0) {
									zz=this[y]-comp[y]
									if(zz>0) {
										return(1) 
									} else {
										if(zz<0) return(-1) else return(0)
									}
								 } else {
									return(0)
								 }
							}
						)
						x.res=mean(x.res[x.res!=0])
						if(is.na(x.res)) x.res=0
					} else {
						yy=this-comp
						zz=yy
						zz[yy>0]=1
						zz[yy<0]=-1
						zz[yy==0]=0
						x.res=mean(zz[zz!=0])
						if(is.na(x.res)) x.res=0
					}
					if(x%in%names(shifts.res$branch.shift.probability)) return(x.res) else return(0)
				}
		)
		colors.nodes=ccx[match(round(shift.direction,1),c.seq)]
		names(colors.nodes)=all.nodes
		colors.nodes=colors.nodes[match(phy$edge[,2], as.numeric(names(colors.nodes)) )]
	} else {
		colors.nodes=NULL
		shift.direction=rep(0,nrow(phy$edge))
	}
	
	## PLOTTING OF TREE ##
	CEX=max(c(0.05, (2+-0.36*log(Ntip(phy))) ))	
	
	if(legend) {
		def.par <- par(no.readonly = TRUE) 
		on.exit(par(def.par))
		layout(matrix(c(1,2,1,3,1,4), 3, 2, byrow=TRUE), widths=c(20,5), respect=FALSE)
	}
	
	if(pdf) pdf(file=paste(outdir, paste(out.base,"pdf",sep="."),sep="/"))
	plot(phy, cex=CEX, lwd=0.4, edge.color=colors.branches, show.tip.label=TRUE, label.offset=strwidth(par("pch"),cex=1.25), no.margin=TRUE, ...)
	NN=phy$edge[,2]
	ll<-cc<-rr<-rep(0,length(NN))
	if(!is.null(shifts.res$branch.shift.probability)) {
		branches=as.numeric(names(shifts.res$branch.shift.probability))
		ll[match(branches, NN)]=1
		cc[match(branches, NN)]=shifts.res$branch.shift.probability
		rr=colors.nodes
	}
	
	edgelabels.rj(text=NULL, pch=ifelse(ll==1, 21, NA), cex=4*cc, col=add.transparency(rr,0.95), bg=add.transparency(rr,0.5), lwd=0.5)
	## END PLOTTING of TREE ##
	
	if(legend) {
		legend.seq=seq(1,color.length,by=2)
		point.seq=rev(seq(0.2,0.8,length=length(legend.seq)))
		
		# shift direction
		plot(rep(-0.5, length(point.seq)), point.seq, xlim=c(-1,2), ylim=c(0,1), cex=2, pch=21, col = add.transparency(rev(ccx[legend.seq]),0.95), bg = add.transparency(rev(ccx[legend.seq]),0.5), bty="n", xaxt="n", yaxt="n")
		mtext("shift direction",side=3,line=-3,cex=0.75)
		text(rep(1, length(point.seq)), point.seq, adj=1, labels=sprintf("%9.2f", rev(c.seq[legend.seq])))

		# posterior estimates
		if(any(colors.branches!=1)) {
			cbt=colors.branches.tmp$legend.seq
			lchars=sapply(cbt, floor)
			if(all(lchars>0)) {
				ldec=1
			} else if(all(lchars==0)) {
				ldec=min(5, max(nchar(range(cbt))))
			} else {
				ldec=3
			}
				
			lnchar=max(nchar(lchars))+ldec
			plot(rep(-0.5, length(point.seq)), point.seq, xlim=c(-1,2), ylim=c(0,1), cex=2, pch=22, col = "darkgray", bg = rev(colors.branches.tmp$legend.col[legend.seq]), bty="n", xaxt="n", yaxt="n")
			mtext(paste("posterior ",param,sep=""),side=3,line=-3,cex=0.75)
			text(rep(1, length(point.seq)), point.seq, adj=1, labels=sprintf(paste("%",max(10,lnchar),".",ldec,"f",sep=""), cbt[legend.seq]))

		}

		# posterior probabilities of shift
		plot(rep(-0.5, length(point.seq)), point.seq, xlim=c(-1,2), ylim=c(0,1), cex=4*(seq(1, 0, length=9)), pch=21, col = add.transparency("darkgray",0.8), bg = add.transparency("white",0.8), bty="n", xaxt="n", yaxt="n")
		mtext("shift probability",side=3,line=-3,cex=0.75)
		text(rep(1, length(point.seq)), point.seq, adj=1, labels=sprintf("%10.3f", seq(1, 0, length=9)))
		
		# reset plotting device
		invisible()
	}
	if(pdf) dev.off()
	
	# GENERATE TABULAR OUTPUT of RESULTS
	if(verbose) {
		if(length(!is.null(shifts.res$descendants))) {
		    dd=as.numeric(names(shifts.res$descendants))
			for(nn in 1:length(dd)) {
				write.table(c(dd[nn], shifts.res$descendants[[nn]]), file=paste(outdir, paste(lab,"shift.descendants.txt",sep="."), sep="/"), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
			}
		}
		
		allres=data.frame(matrix(NA, nrow=nrow(phy$edge), ncol=3))
		shift.direction=shift.direction[match(phy$edge[,2],all.nodes)]
		shift.probability=cc
		allres=data.frame(branch=phy$edge[,2],shift.direction,shift.probability)
		if(nrow(allres)>0) {
			cat("\n\n", rep(" ", 10), toupper("posterior summary"), "\n")
			cat("\n")
			table.print(allres, c(0,4,4), 0)
			names(allres)=c("branch", "relative.shift.direction", "conditional.shift.probability")
			write.table(allres, file=paste(outdir, paste(lab, "estimates.summary.txt", sep="."), sep="/"), quote=FALSE, row.names=FALSE)
		}
	}
	setwd(oldwd)
	return(list(estimates=ests, shifts=shifts))
}

#general phylogenetic plotting utility, given a named vector or data.frame of values that can be associated with phy$edge[,2]
#author: JM EASTMAN 2010
#note: small values are given bluish hues, large values reddish hues; median values are given gray hues

branchcol.plot <-
function(phy, cur.rates, color.length=8, digits=3, plot=TRUE, legend=TRUE, legend.title="", log=FALSE, ...) {
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



#rjmcmc run diagnosis, generating tallies of proposed and accepted updates by class of proposal mechanism
#author: JM EASTMAN 2010

summarize.run <-
function(n.accept, n.props, prop.names) {
	df=data.frame(cbind(proposed=n.props, accepted=n.accept, adoptrate=n.accept/n.props))
	rownames(df)=prop.names
	cat("\n\n",rep(" ",10),toupper(" sampling summary"),"\n")
	
	if(any(is.na(df))) df[is.na(df)]=0
	table.print(df, digits=c(0,0,4), buffer=6)
	
	cat("\n\n")
	
}

#general phylogenetic plotting utility, which is a modification of ape:::edgelabels, plotting edge symbols at the temporal beginning of the branch
#author: E PARADIS 2009 and JM EASTMAN 2010
#note: may not be trustworthy where lastPP$type is not "phylogram"

edgelabels.rj <-
function (text, edge, adj = c(0.5, 0.5), frame = "rect", pch = NULL, 
thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
bg = "lightgreen", horiz = FALSE, width = NULL, height = NULL, 
...) 
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(edge)) {
        sel <- 1:dim(lastPP$edge)[1]
        subedge <- lastPP$edge
    }
    else {
        sel <- edge
        subedge <- lastPP$edge[sel, , drop = FALSE]
    }
    if (lastPP$type == "phylogram") {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            XX <- (lastPP$xx[subedge[, 1]])
            YY <- lastPP$yy[subedge[, 2]]
        }
        else {
            XX <- lastPP$xx[subedge[, 2]]
            YY <- (lastPP$yy[subedge[, 1]])
        }
    }
    else {
        XX <- (lastPP$xx[subedge[, 2]])
        YY <- (lastPP$yy[subedge[, 2]])
    }
	if(missing(text)) text=lastPP$edge[,2]
	
    ape:::BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo, pie, 
					 piecol, col, bg, horiz, width, height, ...)
}



#general printing utility for ensuring equal numbers of characters within columns and defining spacing between columns
#author: JM EASTMAN 2010
#note: works only for numeric dataframes

table.print=function(df,digits=4,buffer=5){
	if(length(buffer) != ncol(df) | length(buffer)==1) buffer=rep(buffer[1],ncol(df))
	if(length(digits) != ncol(df) | length(digits)==1) digits=rep(digits[1],ncol(df))
	ss=sapply(round(df),nchar)
	lar=df>1
	nn=sapply(names(df),nchar)
	
# find longest string
	strw=sapply(1:ncol(df), function(x) max(nn, max(1,(ss[lar])+digits[x],na.rm=TRUE),na.rm=TRUE))  
	pr.df=data.frame(sapply(1:ncol(df), function(x) sprintf(paste("%",(strw[x]+buffer[x]),".",digits[x],"f",sep=""),df[,x])))   
	names(pr.df)=names(df)
	rownames(pr.df)=rownames(df)
	print(pr.df)
}


#generic plotting utility for adding transparency to a RGB color (specified by a 7 or 9 character string)
#author: R FITZ-JOHN 2010

add.transparency <-
function (col, alpha) {
    tmp <- col2rgb(col)/255
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
}



