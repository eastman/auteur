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



#logging utility used by rjmcmc
#author: JM EASTMAN 2010

parlogger <- function(phy, init=FALSE, primary.parameter, parameters, parmBase, runLog, end=FALSE, class=c("rj","hmm")) {
	if(missing(class)) class="rj"
	
	# determine if rjMCMC or hmmMCMC: parlogs are edgewise parameters
	if(class=="hmm") {
		parlogs=paste(parmBase, paste(c(primary.parameter), "txt", sep="."), sep="")
	} else {
		parlogs=paste(parmBase, paste(c("shifts",primary.parameter), "txt", sep="."), sep="")
	}
	
	# add eol to files if terminating run
	if(end==TRUE) {
		sapply(1:length(parlogs), function(x) write.table("\n", parlogs[x], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE))
		return()
	}
	
	parameters->parms
	parnames=names(parms)
	ii=match(c("gen"),parnames)
	general=match(c("lnL","lnL.p","lnL.h"),parnames)
	target=match(primary.parameter, parnames)

	delta=match("delta",parnames)
	accessory=match(parnames[-c(ii,delta,target,general)],parnames)
	
	if(init) {
		if(class=="hmm") {
			sapply(parlogs[1:length(parlogs)], function(x) write.table("", file=x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			write.table(paste("state", primary.parameter, paste(parnames[accessory], collapse="\t"), paste(parnames[general],collapse="\t"), sep="\t"), file=runLog, quote=FALSE, col.names=FALSE, row.names=FALSE)		
		} else {
			sapply(parlogs[1:length(parlogs)], function(x) write.table("", file=x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			write.table(paste("state", "min", "max", "mean", primary.parameter, paste(parnames[accessory],collapse="\t"), paste(parnames[general],collapse="\t"), sep="\t"), file=runLog, quote=FALSE, col.names=FALSE, row.names=FALSE)		
		}
	} else {
		dd=parms[[delta]]

		if(class=="hmm") {
			# to log file
			msg<-paste(parms[[ii]], parms[[target]], paste(sprintf("%.3f", parms[accessory]),collapse="\t"), paste(sprintf("%.2f", parms[general]),collapse="\t"),sep="\t")
			write(msg, file=runLog, append=TRUE)
			
			# to parameter files
			out=list(dd)		
			sapply(1:length(parlogs), function(x) write.table(paste(out[[x]],collapse=" "), parlogs[x], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE))
			
		} else {
			pp=parms[[target]]
			K=sum(dd)+1
			range.p=range(pp)
			mean.p=mean(pp)
			
			# to log file
			msg<-paste(parms[[ii]], sprintf("%.3f", range.p[1]), sprintf("%.3f", range.p[2]), sprintf("%.3f", mean.p), K, paste(sprintf("%.3f", parms[accessory]),collapse="\t"), paste(sprintf("%.2f", parms[general]),collapse="\t"),sep="\t")
			write(msg, file=runLog, append=TRUE)
			
			# to parameter files
			out=list(dd,pp)		
			sapply(1:length(parlogs), function(x) write.table(paste(out[[x]],collapse=" "), parlogs[x], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE))
			
		}
	}
}


#rjmcmc utility for initiating a proposal width for Markov sampling
#author: JM EASTMAN 2010
#modified: 02.26.2011 to use spline() in computation of best prop.width

calibrate.proposalwidth <-
function(phy, dat, nsteps=100, model=c("BM","jumpBM","OU"), widths=NULL) {
	if(!withinrange(nsteps, 100, 1000)) {
		if(nsteps>1000) nsteps=1000 else nsteps=100
	}
	if(is.null(widths)) widths=2^(-2:4)
	if(model=="BM") {
		acceptance.rates=sapply(widths, function(x) rjmcmc.bm(phy=phy, dat=dat, ngen=nsteps, prop.width=x, summary=FALSE, fileBase="propwidth")$acceptance.rate)
	} else if(model=="jumpBM") {
#		acceptance.rates=sapply(widths, function(x) mcmc.levy(phy=phy, dat=dat, ngen=nsteps, prop.width=x, summary=FALSE, fileBase="propwidth")$acceptance.rate)
	} else if(model=="OU") {
#		while(1){
#			acceptance.rates=sapply(widths, function(x) {
#								f=try(x<-rjmcmc.ou(phy=phy, dat=dat, ngen=nsteps, prop.width=x, summary=FALSE, fileBase="propwidth")$acceptance.rate,silent=TRUE)
#								if(inherits(f,"try-error")) return(0) else return(x)
#								}
#								)
#			if(any(acceptance.rates!=0))break()
#		}
	}
	s=spline(log(widths,base=2),acceptance.rates)
	m.s=mean(s$y)
	s.pick=which(abs(s$y-m.s)==min(abs(s$y-m.s)))
	return(2^mean(s$x[s.pick]))
}



#utility for converting text files to .rda to compress output from rjmcmc
#author: JM EASTMAN 2011

cleanup.files <-
function(parmBase, model, fileBase, phy){
	if(model=="BM") {
		parm="rates"
		parms=c(parm,"shifts")
	} else if(model=="OU") {
		parm="optima"
		parms=c(parm,"shifts")
	} else if(model=="jumpBM") {
		parms="jumps"
	}
	
# read in raw data
	posteriorsamples=lapply(parms, function(x) {
			  y=read.table(paste(parmBase,paste(x,"txt",sep="."),sep="/"))
			  names(y)=phy$edge[,2]
			  return(y)
			  })
	names(posteriorsamples)=parms
	
	save(posteriorsamples, file=paste(parmBase,paste(fileBase,"posteriorsamples","rda",sep="."),sep="/"))
	unlink(paste(parmBase,paste(parms,"txt",sep="."),sep="/"))
}


#general phylogenetic statistical utility for comparing values ('posterior.rates') from posterior densities for two sets of supplied branches ('nodes')
#author: JM EASTMAN 2010

compare.rates <-
function(branches=list(A=c(NULL,NULL),B=c(NULL,NULL)), phy, posterior.values, ...){
	nodes=branches
	rates=posterior.values
	if(is.null(names(nodes))) names(nodes)=c("A","B")
	rates.a=rates[,match(nodes[[1]],names(rates))]
	rates.b=rates[,match(nodes[[2]],names(rates))]
	edges.a=phy$edge.length[match(nodes[[1]],phy$edge[,2])]
	edges.b=phy$edge.length[match(nodes[[2]],phy$edge[,2])]
	weights.a=edges.a/sum(edges.a)
	weights.b=edges.b/sum(edges.b)
	if(length(nodes[[1]])>1) {r.scl.a=apply(rates.a, 1,function(x) sum(x*weights.a))} else {r.scl.a=sapply(rates.a, function(x) sum(x*weights.a))}
	if(length(nodes[[2]])>1) r.scl.b=apply(rates.b, 1,function(x) sum(x*weights.b)) else r.scl.b=sapply(rates.b, function(x) sum(x*weights.b))
	return(list(scl.A=r.scl.a, scl.B=r.scl.b, r.test=randomization.test(r.scl.a, r.scl.b, ...)))	
}

#general phylogenetic utility for rapid computation of the maximum-likelihood estimate of the Brownian motion rate parameter (unless there are polytomies, in which case the function wraps geiger:::fitContinuous)
#author: L REVELL 2010, LJ HARMON 2009, and JM EASTMAN 2010

fit.continuous <-
function(phy, dat){
	n=Ntip(phy)
	ic=try(pic(dat,phy),silent=TRUE)
	if(!inherits(ic, "try-error")) {
		r=mean(ic^2)*((n-1)/n)
		return(r)
	} else {
		return(fitContinuous(phy,dat)$Trait1$beta)
	}
}

#general phylogenetic utility to find starting point of alpha and sigmasq under an OU process
#author: JM EASTMAN 2011

fit.hansen <-
function(phy, dat, alphabounds){
	ot=ape2ouch(phy)

	otd <- as(ot,"data.frame")
	ot <- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
	otd$regimes <- as.factor(1)
	dat.df<-data.frame(trait=dat,row.names=names(dat))
	dat.df$labels <- rownames(dat.df)
	otd <- merge(otd,dat.df,by="labels",all=TRUE)
	rownames(otd) <- otd$nodes
	h=hansen(tree=ot, data=otd[c("trait")], regimes=otd["regimes"],sqrt.alpha=1,sigma=1)
	alpha=h@sqrt.alpha^2 
	sigmasq=h@sigma^2
	if(withinrange(alpha, min(alphabounds), max(alphabounds))) alpha=alpha else alpha=alphabounds[which(abs(alphabounds-alpha)==min(abs(alphabounds-alpha)))]
	return(list(alpha=h@sqrt.alpha^2, sigmasq=h@sigma^2))
}


#logging utility used by rjmcmc
#author: JM EASTMAN 2010

generate.error.message <-
function(ntip, i, mod.cur, mod.new, lnL, lnPrior, lnHastings, curCats, newCats, errorLog) {
	if(!file.exists(file=errorLog)) {
		write(paste("gen", "cur.lnL", "new.lnL", "lnL", "lnPrior", "lnHastings", "cur.catg", "new.catg", sep="\t"), file=errorLog)
	} 
	write(paste(i, sprintf("%.3f", mod.cur$lnL), sprintf("%.3f", mod.new$lnL), sprintf("%.3f", lnL), sprintf("%.3f", lnPrior), sprintf("%.3f", lnHastings), sum(curCats), sum(newCats), sep="\t"),  file=errorLog, append=TRUE)
	
}

#utility for generating a list of indicator variables and values for each branch in a phylogeny, using a random model complexity as a starting point (or by constraining model complexity)
#author: JM EASTMAN 2010

generate.starting.point <-
function(data, phy, node.des=NULL, logspace=TRUE, K=FALSE, prop.width) { 

	if(is.null(node.des)) {
		node.des		<- sapply(unique(c(phy$edge[1,1],phy$edge[,2])), function(x) get.descendants.of.node(x, phy))
		names(node.des) <- c(phy$edge[1,1], unique(phy$edge[,2]))
	}
	
	nn=length(phy$edge.length)
	ntip<-n<-Ntip(phy)
	if(!K) nshifts=rtpois(1,log(ntip),nn) else nshifts=K-1
	if(nshifts!=0) bb=sample(c(rep(1,nshifts),rep(0,nn-nshifts)),replace=FALSE) else bb=rep(0,nn)
	names(bb)=phy$edge[,2]
	shifts=as.numeric(names(bb)[bb==1])
	if(!logspace) {
		init.rate=mean(data)
		vd=sapply(shifts, function(x) {
				  z=c()
				  xx=unlist(node.des[[which(names(node.des)==x)]])
				  xx=xx[xx<=n]
				  z=c(xx,x[x<=n])
				  yy=c(mean(data[z]))
				return(yy)})
		values=c(vd, init.rate)
	} else {
		init.rate=fit.continuous(phy,data)
		min.max=c(adjustrate(init.rate, prop.width), adjustrate(init.rate, prop.width))
		values=runif(sum(bb)+1, min=min.max[min.max==min(min.max)], max=min.max[min.max==max(min.max)])
	}
	internal.shifts<-tip.shifts<-numeric(0)
	internal.shifts=sort(shifts[shifts>ntip])
	tip.shifts=shifts[shifts<=ntip]
	
	if(length(internal.shifts)==0 & length(tip.shifts)==0) {
		vv=rep(values, nn)
	} else {
		vv=bb
		vv[]=values[length(values)]
		i=0
		if(length(internal.shifts)!=0) {
			for(i in 1:length(internal.shifts)) {
				d=node.des[which(names(node.des)==internal.shifts[i])]
				vv[match(c(internal.shifts[i], unlist(d)), names(vv))]=values[i]
			}
		}
		if(length(tip.shifts)!=0) {
			for(j in 1:length(tip.shifts)) {
				vv[match(tip.shifts[j], names(vv))]=values[j+i]
			}
		}
	}
	return(list(delta=unname(bb), values=unname(vv)))
}

#general utility for combining supplied list of dataframes into a single 'intercalated' sample. E.g., list(as.data.frame(c(AAA)),as.data.frame(c(BBB))) becomes data.frame(c(ABABAB))
#author: JM EASTMAN 2010

intercalate.samples <-
function(list.obj) {
	if(length(list.obj)<2) warning("Nothing to be done... too few objects in list.")
	vv=sapply(list.obj, is.vector)
	if(all(vv)) list.obj=lapply(list.obj, function(x) {y=data.frame(x); names(y)=NULL; return(y)})
	orig.names=names(list.obj[[1]])
	n=length(list.obj)
	R=sapply(list.obj, nrow)
	if(length(unique(R))!=1) {
		minR=nrow(list.obj[[which(R==min(R))]])
		list.obj=lapply(list.obj, function(x) {y=as.data.frame(x[1:minR,]); names(y)=orig.names; return(y)})
		warning("Objects pruned to the length of the smallest.")
	}
	C=sapply(list.obj, ncol)
	if(length(unique(C))!=1) stop("Cannot process objects; column lengths do not match.")
	c=unique(C)
	if(c>1) {
		M=lapply(1:length(list.obj), function(x) mm=match(names(list.obj[[x]]), orig.names))
		for(mm in 2:length(list.obj)) {
			list.obj[[mm]]=list.obj[[mm]][,M[[mm]]]
		}
	}
	r=nrow(list.obj[[1]])
	indices=lapply(1:n, function(x) seq(x, r*n, by=n))
	
	out.array=array(dim=c(r*n, c))
	for(nn in 1:n){
		out.array[indices[[nn]],]=as.matrix(list.obj[[nn]])
	}
	
	out.array=data.frame(out.array)
	if(!is.null(orig.names)) names(out.array)=orig.names else names(out.array)=paste("X",1:ncol(out.array),sep="")
	return(out.array)
}



#intercalates samples from multiple independent Markov chains (generated by rjmcmc.bm())
#author: JM EASTMAN 2010

pool.rjmcmcsamples <-
function(base.dirs, lab="combined"){
	outdir=paste(lab,"combined.rjmcmc",sep=".")
	if(!file.exists(outdir)) dir.create(outdir)
	
# estimated parameters
	samples=lapply(base.dirs, function(x) {load(paste(x, dir(x, pattern="posteriorsamples.rda"), sep="/")); return(posteriorsamples)})
	nn=unique(unlist(lapply(samples, names)))
	posteriorsamples=list()
	for(i in nn) {
		res.sub=intercalate.samples(lapply(1:length(samples), function(x) samples[[x]][[i]]))
		posteriorsamples[[which(nn==i)]]=res.sub
	}
	names(posteriorsamples)=nn
	save(posteriorsamples, file=paste(outdir, paste(lab,"posteriorsamples.rda",sep="."), sep="/"))
	
# logfile
	
	write.table(intercalate.samples(lapply(base.dirs, function(x) {y=read.table(paste(x, dir(x, pattern="rjmcmc.log"), sep="/"), header=TRUE); return(y)})), paste(outdir, paste(lab,"rjmcmc.log",sep="."), sep="/"), quote=FALSE, row.names=FALSE, sep="\t")
}

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



