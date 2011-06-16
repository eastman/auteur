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
	if(length(nn) & any(!is.na(nn))) {
		desc=lapply(nn, function(x) {y=get.descendants.of.node(x, phy=phy, tips=TRUE); if(length(y)) return(phy$tip.label[y]) else return(NULL)})
		names(desc)=nn		
		di=sapply(nn, function(x) if(x>Ntip(phy)) return(TRUE) else return(FALSE))
		if(any(di)) desc=desc[di] else desc=NULL
	} else {
		desc<-branches<-NULL
	}
	return(list(descendants=desc, branch.shift.probability=branches))
}

shifts.plot <-
function(phy, base.dir, burnin=0, level=0.01, internal.only=FALSE, paint.branches=TRUE, legend=TRUE, verbose=TRUE, lab="OUT", ...) {
	color.length=17
	oldwd=getwd()
	setwd(base.dir)
	
	if(verbose) {
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
	
#	if(pdf) pdf(file=paste(outdir, paste(out.base,"pdf",sep="."),sep="/"))
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
#	if(pdf) dev.off()
	
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

## from ape:::plot.phylogram
plotfangram=function (x, type = "fan", use.edge.length = TRUE, node.pos = NULL, 
show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
direction = "rightwards", lab4ut = "horizontal", tip.color = "black", default=list(col="black", lwd=1), 
...) 
{
    Ntip <- length(x$tip.label)
    if (Ntip == 1) {
        warning("found only one tip in the tree")
        return(NULL)
    }
    if (any(tabulate(x$edge[, 1]) == 1)) 
	stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C("node_height", 
															 as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
																											 1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
															 DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C("node_depth", 
														as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
																										1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
																																							  Nnode), DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
									 edge.length) .C("node_depth_edgelength", as.integer(Ntip), 
													 as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
																										  2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
																																								 Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
							  "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
										"upwards", "downwards"))
    if (is.null(x$edge.length)) 
	use.edge.length <- FALSE
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
	root.edge <- FALSE
    if (type == "fan" && root.edge) {
        warning("drawing root edge with type = 'fan' is not yet supported")
        root.edge <- FALSE
    }
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge
    if (phyloORclado) {
        phyOrder <- attr(x, "order")
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x)
            if (!identical(x$edge, xe)) {
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
					edge.color <- rep(edge.color, length.out = Nedge)
					edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
					edge.width <- rep(edge.width, length.out = Nedge)
					edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
					edge.lty <- rep(edge.lty, length.out = Nedge)
					edge.lty <- edge.lty[ereorder]
                }
            }
        }
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = "pruningwise")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) 
			node.pos <- 2
        }
        if (node.pos == 1) 
		yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
        else {
            ans <- .C("node_height_clado", as.integer(Ntip), 
					  as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 
																			 2]), as.integer(Nedge), double(Ntip + Nnode), 
					  as.double(yy), DUP = FALSE, PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
			xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
			1
            xx <- max(xx) - xx
        }
        else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
									   z$edge.length)
        }
    }
    else switch(type, fan = {
				TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
				xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
				theta <- double(Ntip)
				theta[TIPS] <- xx
				theta <- c(theta, numeric(Nnode))
				theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, theta)
				if (use.edge.length) {
				r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
										  z$edge.length)
				} else {
				r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
				r <- 1/r
				}
				xx <- r * cos(theta)
				yy <- r * sin(theta)
				}, unrooted = {
				nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
				XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, z$edge, 
													   z$edge.length, nb.sp) else unrooted.xy(Ntip, Nnode, 
																							  z$edge, rep(1, Nedge), nb.sp)
				xx <- XY$M[, 1] - min(XY$M[, 1])
				yy <- XY$M[, 2] - min(XY$M[, 2])
				}, radial = {
				X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
				X[X == 1] <- 0
				X <- 1 - X/Ntip
				yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
				Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
				xx <- X * cos(Y)
				yy <- X * sin(Y)
				})
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") 
			xx <- xx + x$root.edge
            if (direction == "upwards") 
			yy <- yy + x$root.edge
        }
    }
    if (no.margin) 
	par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par("pin")[1]
                strWi <- strwidth(x$tip.label, "inches")
                xx.tips <- xx[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * xx.tips + 
												   strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
				tmp <- max(xx.tips) * 1.5
                else {
					tmp <- if (show.tip.label) 
                    max(xx.tips + strWi/alp)
					else max(xx.tips)
                }
                x.lim[2] <- tmp
            }
            else x.lim <- c(1, Ntip)
        }
        else switch(type, fan = {
					if (show.tip.label) {
					offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
								  cex)
					x.lim <- c(min(xx) - offset, max(xx) + offset)
					} else x.lim <- c(min(xx), max(xx))
					}, unrooted = {
					if (show.tip.label) {
					offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
								  cex)
					x.lim <- c(0 - offset, max(xx) + offset)
					} else x.lim <- c(0, max(xx))
					}, radial = {
					if (show.tip.label) {
					offset <- max(nchar(x$tip.label) * 0.03 * cex)
					x.lim <- c(-1 - offset, 1 + offset)
					} else x.lim <- c(-1, 1)
					})
    }
    else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
		x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
		x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
						 cex)
        if (type == "radial") 
		x.lim[1] <- if (show.tip.label) 
		-1 - max(nchar(x$tip.label) * 0.03 * cex)
		else -1
    }
    if (phyloORclado && direction == "leftwards") 
	xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
			y.lim <- c(1, Ntip)
            else {
                y.lim <- c(0, NA)
                pin2 <- par("pin")[2]
                strWi <- strwidth(x$tip.label, "inches")
                yy.tips <- yy[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * yy.tips + 
												   strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
				tmp <- max(yy.tips) * 1.5
                else {
					tmp <- if (show.tip.label) 
                    max(yy.tips + strWi/alp)
					else max(yy.tips)
                }
                y.lim[2] <- tmp
            }
        }
        else switch(type, fan = {
					if (show.tip.label) {
					offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
								  cex)
					y.lim <- c(min(yy) - offset, max(yy) + offset)
					} else y.lim <- c(min(yy), max(yy))
					}, unrooted = {
					if (show.tip.label) {
					offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
								  cex)
					y.lim <- c(0 - offset, max(yy) + offset)
					} else y.lim <- c(0, max(yy))
					}, radial = {
					if (show.tip.label) {
					offset <- max(nchar(x$tip.label) * 0.03 * cex)
					y.lim <- c(-1 - offset, 1 + offset)
					} else y.lim <- c(-1, 1)
					})
    }
    else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
		y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
		y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
						 cex)
        if (type == "radial") 
		y.lim[1] <- if (show.tip.label) 
		-1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
				 cex)
		else -1
    }
    if (phyloORclado && direction == "downwards") 
	yy <- y.lim[2] - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
		x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
		y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 
	1
    else NA
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, 
		 axes = FALSE, asp = asp, ...)
    if (is.null(adj)) 
	adj <- if (phyloORclado && direction == "leftwards") 
	1
	else 0
    if (phyloORclado && show.tip.label) {
        MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
        loy <- 0
        if (direction == "rightwards") {
            lox <- label.offset + MAXSTRING * 1.05 * adj
        }
        if (direction == "leftwards") {
            lox <- -label.offset - MAXSTRING * 1.05 * (1 - adj)
        }
        if (!horizontal) {
            psr <- par("usr")
            MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
															   psr[1])
            loy <- label.offset + MAXSTRING * 1.05 * adj
            lox <- 0
            srt <- 90 + srt
            if (direction == "downwards") {
                loy <- -loy
                srt <- 180 + srt
            }
        }
    }
    if (type == "phylogram") {
        phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
					   edge.color, edge.width, edge.lty)
    }
    else {
        if (type == "fan") {
            ereorder <- match(z$edge[, 2], x$edge[, 2])
            if (length(edge.color) > 1) {
                edge.color <- rep(edge.color, length.out = Nedge)
                edge.color <- edge.color[ereorder]
            }
            if (length(edge.width) > 1) {
                edge.width <- rep(edge.width, length.out = Nedge)
                edge.width <- edge.width[ereorder]
            }
            if (length(edge.lty) > 1) {
                edge.lty <- rep(edge.lty, length.out = Nedge)
                edge.lty <- edge.lty[ereorder]
            }
            circularplot(z$edge, Ntip, Nnode, xx, yy, theta, 
						 r, edge.color, edge.width, edge.lty, default=default) #
        }
        else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
							edge.lty)
    }
    if (root.edge) 
	switch(direction, rightwards = segments(0, yy[ROOT], 
											x$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
																						 yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]), upwards = segments(xx[ROOT], 
																																						 0, xx[ROOT], x$root.edge), downwards = segments(xx[ROOT], 
																																																		 yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge))
    if (show.tip.label) {
        if (is.expression(x$tip.label)) 
		underscore <- TRUE
        if (!underscore) 
		x$tip.label <- gsub("_", " ", x$tip.label)
        if (phyloORclado) 
		text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
			 adj = adj, font = font, srt = srt, cex = cex, 
			 col = tip.color)
        if (type == "unrooted") {
            if (lab4ut == "horizontal") {
                y.adj <- x.adj <- numeric(Ntip)
                sel <- abs(XY$axe) > 0.75 * pi
                x.adj[sel] <- -strwidth(x$tip.label)[sel] * 1.05
                sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
				pi
                x.adj[sel] <- -strwidth(x$tip.label)[sel] * (2 * 
															 abs(XY$axe)[sel]/pi - 0.5)
                sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
                y.adj[sel] <- strheight(x$tip.label)[sel]/2
                sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
                y.adj[sel] <- -strheight(x$tip.label)[sel] * 
				0.75
                text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + y.adj * 
					 cex, x$tip.label, adj = c(adj, 0), font = font, 
					 srt = srt, cex = cex, col = tip.color)
            }
            else {
                adj <- abs(XY$axe) > pi/2
                srt <- 180 * XY$axe/pi
                srt[adj] <- srt[adj] - 180
                adj <- as.numeric(adj)
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                if (label.offset) {
					xx.tips <- xx.tips + label.offset * cos(XY$axe)
					yy.tips <- yy.tips + label.offset * sin(XY$axe)
                }
                font <- rep(font, length.out = Ntip)
                tip.color <- rep(tip.color, length.out = Ntip)
                cex <- rep(cex, length.out = Ntip)
                for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
									   cex = cex[i], x$tip.label[i], adj = adj[i], 
									   font = font[i], srt = srt[i], col = tip.color[i])
            }
        }
        if (type %in% c("fan", "radial")) {
            xx.tips <- xx[1:Ntip]
            yy.tips <- yy[1:Ntip]
            angle <- atan2(yy.tips, xx.tips)
            if (label.offset) {
                xx.tips <- xx.tips + label.offset * cos(angle)
                yy.tips <- yy.tips + label.offset * sin(angle)
            }
            s <- xx.tips < 0
            angle <- angle * 180/pi
            angle[s] <- angle[s] + 180
            adj <- as.numeric(s)
            font <- rep(font, length.out = Ntip)
            tip.color <- rep(tip.color, length.out = Ntip)
            cex <- rep(cex, length.out = Ntip)
            for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], x$tip.label[i], 
								   font = font[i], cex = cex[i], srt = angle[i], 
								   adj = adj[i], col = tip.color[i])
        }
    }
    if (show.node.label) 
	text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
		 x$node.label, adj = adj, font = font, srt = srt, 
		 cex = cex)
    L <- list(type = type, use.edge.length = use.edge.length, 
			  node.pos = node.pos, show.tip.label = show.tip.label, 
			  show.node.label = show.node.label, font = font, cex = cex, 
			  adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
			  x.lim = x.lim, y.lim = y.lim, direction = direction, 
			  tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
		   envir = .PlotPhyloEnv)
    invisible(L)
}

## FROM ape:::circular.plot
circularplot =function (edge, Ntip, Nnode, xx, yy, theta, r, edge.color, edge.width, 
edge.lty, default=list(col="black", lwd=1)) 
{
    r0 <- r[edge[, 1]]
    r1 <- r[edge[, 2]]
    theta0 <- theta[edge[, 2]]
    costheta0 <- cos(theta0)
    sintheta0 <- sin(theta0)
    x0 <- r0 * costheta0
    y0 <- r0 * sintheta0
    x1 <- r1 * costheta0
    y1 <- r1 * sintheta0
    segments(x0, y0, x1, y1, col = edge.color, lwd = edge.width, 
			 lty = edge.lty)
    tmp <- which(diff(edge[, 1]) != 0)
    start <- c(1, tmp + 1)
    Nedge <- dim(edge)[1]
    end <- c(tmp, Nedge)
    foo <- function(edge.feat, default) {
        if (length(edge.feat) == 1) 
		return(rep(edge.feat, Nnode))
        else {
            edge.feat <- rep(edge.feat, length.out = Nedge)
            feat.arc <- rep(default, Nnode)
            for (k in 1:Nnode) {
                tmp <- edge.feat[start[k]]
                if (tmp == edge.feat[end[k]]) 
				feat.arc[k] <- tmp
            }
        }
        feat.arc
    }
    co <- foo(edge.color, default$col) #
    lw <- foo(edge.width, default$lwd) #
    ly <- foo(edge.lty, 1)
    for (k in 1:Nnode) {
        i <- start[k]
        j <- end[k]
        X <- rep(r[edge[i, 1]], 100)
        Y <- seq(theta[edge[i, 2]], theta[edge[j, 2]], length.out = 100)
        lines(X * cos(Y), X * sin(Y), col = co[k], lwd = lw[k], 
			  lty = ly[k])
    }
}



