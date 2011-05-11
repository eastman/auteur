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
	shifts=posteriorsamples$rate.shifts
	burnin=ceiling(burnin*nrow(shifts))
	shifts=shifts[-c(1:(burnin)),]
	shifts.res=process.shifts(phy, shifts, level)
	ests=posteriorsamples$rates
	ests=ests[-c(1:(burnin)),]
	average.ests.tmp=apply(ests,2,mean)
	average.ests=average.ests.tmp-median(average.ests.tmp)
	
	# collect edge colors (for rates)
	if(paint.branches) {
		colors.branches.tmp=branchcol.plot(phy, ests, plot=FALSE, color.length=color.length)
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

		# posterior rates
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
			mtext("posterior rates",side=3,line=-3,cex=0.75)
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
}

