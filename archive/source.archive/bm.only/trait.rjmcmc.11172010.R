## Written by LJ HARMON (2008), AL HIPP (2008), and JM EASTMAN (2010)

rjMCMC.trait<-function (phy, data, ngen=1000, sampleFreq=100, probMergeSplit=0.05, probRoot=0.01, lambdaK=log(2), constrainK=FALSE, jumpsize=2, fileBase="result", simplestart=FALSE) 
{
	model="BM"
	heat=1
#	require(ouch)
	require(geiger)
		
### prepare data for rjMCMC
	cur.model		<- model								

	dataList		<- prepare.data(phy, data)			
	ape.tre			<- cur.tre <- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat						
	nn				<- length(ape.tre$edge.length)	
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))
	subtrees		<- subtrees(phy)
	subtrees.marker <- sapply(subtrees, function(x) return(min(x$node.label)))

# initialize parameters
	if(is.numeric(constrainK) & (constrainK > length(ape.tre$edge) | constrainK < 1)) stop("Constraint on rate shifts is nonsensical. Ensure that constrainK is at least 1 and less than the number of available nodes in the tree.")

	if(simplestart | is.numeric(constrainK) ) {
		if(is.numeric(constrainK)) {
			init.rate	<- generate.starting.point(data, ape.tre, node.des, theta=FALSE, K=constrainK, jumpsize=jumpsize)
		} else {
			init.rate	<- list(values=rep(fitContinuous(phy,data)$Trait1$beta,length(phy$edge.length)),delta=rep(0,length(phy$edge.length)))
		}
	} else {
		init.rate	<- generate.starting.point(data, ape.tre, node.des, theta=FALSE, K=constrainK, jumpsize=jumpsize )
	}
	
    cur.rates		<- init.rate$values			
	cur.delta.rates	<- init.rate$delta
	cur.root		<- adjust.value(mean(orig.dat), jumpsize)
	cur.tre			<- apply.BMM.to.apetree(ape.tre, cur.rates)
	
	mod.cur = new.bm.lik.fn(cur.rates, cur.root, cur.tre, orig.dat)
	
# proposal counts
	nRateProp =		0	
	nRateSwapProp = 0									
	nRcatProp =		0										
	nRootProp =		0
	nRateOK =		0											
	nRateSwapOK =	0
	nRcatOK =		0											
	nRootOK =		0

	cOK=c(												
		nRateOK,
		nRateSwapOK,
		nRcatOK, 
		nRootOK 
	)
		
	tickerFreq=ceiling(ngen/30)
	
# file handling
	parmBase=paste(model, fileBase, "parameters/",sep=".")
	if(!file.exists(parmBase)) dir.create(parmBase)
	parlogger(model=model, init=TRUE, node.des, parmBase=parmBase)
	errorLog=file(paste(parmBase,paste(cur.model, fileBase, "rjTraitErrors.log",sep="."),sep="/"),open='w+')
	generate.error.message(i=NULL, mod.cur=NULL, mod.new=NULL, lnR=NULL, errorLog=errorLog, init=TRUE)
	runLog=file(paste(parmBase,paste(cur.model, fileBase, "rjTrait.log",sep="."),sep="/"),open='w+')
	generate.log(bundled.parms=NULL, cur.model, file=runLog, init=TRUE)


### Begin rjMCMC
    for (i in 1:ngen) {
        lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
		startparms = c(nRateProp, nRateSwapProp, nRcatProp, nRootProp)

## BM IMPLEMENTATION ##
		while(1) {
			if (runif(1) < (2 * probMergeSplit) & !constrainK) {			# adjust rate categories
				nr=split.or.merge(cur.delta.rates, cur.rates, ape.tre, node.des, lambdaK)
				new.rates=nr$new.values
				new.delta.rates=nr$new.delta
				new.root=cur.root
				nRcatProp=nRcatProp+1
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
				break()
			} else { 
				if(runif(1)<probRoot) {										# adjust root
					new.root=adjust.value(cur.root, jumpsize)
					new.rates=cur.rates
					new.delta.rates=cur.delta.rates
					nRootProp=nRootProp+1
					break()
				} else {													
					if(runif(1)>0.1 & length(unique(cur.rates))>1) {		# swap rates classes
						nr=switcheroo(cur.rates)
						new.rates=nr$new.values ##
						new.delta.rates=cur.delta.rates 
						new.root=cur.root
						nRateSwapProp=nRateSwapProp+1
						break()
					} else {												# adjust rates
						new.rates=adjust.rate(cur.rates, jumpsize)
						new.delta.rates=cur.delta.rates
						new.root=cur.root
						nRateProp=nRateProp+1
						break()
					} 
				}
			}
		}
		
		new.tre=apply.BMM.to.apetree(ape.tre, new.rates)
		mod.new=try(new.bm.lik.fn(new.rates, new.root, new.tre, orig.dat),silent=TRUE)
		
		if(inherits(mod.new, "try-error")) {mod.new=as.list(mod.new); mod.new$lnL=Inf}
		if(!is.infinite(mod.new$lnL)) {
			lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
		} else {
			lnLikelihoodRatio = -Inf
			failures=paste("./failed", model, "parameters/",sep=".")
			failed.trees=paste(failures, "trees", sep="/")
			failed.plots=paste(failures, "pdf", sep="/")
			lapply(list(failures, failed.trees, failed.plots), function(x) if(!file.exists(x)) dir.create(x))
				
			pdf(file=paste(failed.plots, paste(i,fileBase,"failed.model.pdf",sep="."),sep="/"))
				plot(new.tre, tip=FALSE)
			dev.off()
			
			write.tree(new.tre, file=paste(failed.trees, paste(i,fileBase,"tre",sep="."),sep="/"))
			write.table(matrix(c(i, new.rates),nrow=1),file=paste(failures,"failed.rates.txt",sep="/"),append=TRUE,col=FALSE,row=FALSE,quote=FALSE)
			write.table(matrix(c(i, new.root),nrow=1),file=paste(failures,"failed.root.txt",sep="/"),append=TRUE,col=FALSE,row=FALSE,quote=FALSE)
		}
		
	# compare likelihoods
		endparms = c(nRateProp=nRateProp, nRateSwapProp=nRateSwapProp, nRcatProp=nRcatProp, nRootProp=nRootProp)

		r=assess.lnR((heat * lnLikelihoodRatio + heat * lnPriorRatio + lnHastingsRatio)->lnR)

		if (runif(1) <= r$r) {			# adopt proposal
			decision="adopt"
			cur.root <- new.root
			cur.rates <- new.rates
			cur.delta.rates <- new.delta.rates
			curr.lnL <- mod.new$lnL
			cur.tre <- new.tre
			mod.cur <- mod.new
			cOK <- determine.accepted.proposal(startparms, endparms, cOK)
		} else {						# deny proposal	
			decision="reject"
			curr.lnL <- mod.cur$lnL 
		}
		
		if(is.infinite(curr.lnL)) stop("starting point has exceptionally poor likelihood")

		if(r$error) generate.error.message(i, mod.cur, mod.new, lnR, errorLog)
		
		if(i%%tickerFreq==0) {
			if(i==tickerFreq) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
			cat(". ")
		}
		
		if(i%%sampleFreq==0) {
			bundled.parms=list(gen=i, mrate=exp(mean(log(cur.rates))), cats=sum(cur.delta.rates)+1, root=cur.root, lnL=curr.lnL)
			generate.log(bundled.parms, cur.model, file=runLog)			
			parlogger(model=model, init=FALSE, node.des=node.des, i=i, curr.lnL=curr.lnL, cur.root=cur.root, cur.rates=cur.rates, cur.delta.rates=cur.delta.rates, parmBase=parmBase)
		}
	}
	
# End rjMCMC
	close(errorLog)
	close(runLog)
	summarize.run(cOK, endparms, cur.model)	
}


## AUXILIARY FUNCTIONS

new.bm.lik.fn <- function(rates, root, ape.tre, orig.dat) { # from LJ HARMON	
# mod 10.18.2010 JM Eastman: using determinant()$modulus rather than det() to stay in log-space
	tree=ape.tre
	y=orig.dat
	
	b <- vcv.phylo(ape.tre)
	w <- rep(root, nrow(b))
	num <- -t(y-w)%*%solve(b)%*%(y-w)/2
	den <- 0.5*(length(y)*log(2*pi) + as.numeric(determinant(b)$modulus))
	list(
		 root=root,
		 lnL = (num-den)[1,1]
		 )	
}

split.or.merge <- function(cur.delta, cur.values, phy, node.des, lambda=lambda) { 
# mod 11.17.2010 JM Eastman
#	cur.delta=cur.delta.rates; cur.values=cur.rates; theta=FALSE; lambda=log(2); jumpsize=2
	bb=cur.delta
	vv=cur.values
	names(vv)<-names(bb)<-phy$edge[,2]
	new.bb=choose.one(bb)
	new.vv=vv
	
	s=names(new.bb[bb!=new.bb])
	all.shifts=as.numeric(names(bb[bb>0]))
	all.D=node.des[[which(names(node.des)==s)]]
	if(length(all.D)!=0) {
		untouchable=unlist(lapply(all.shifts[all.shifts>s], function(x)node.des[[which(names(node.des)==x)]]))
		remain.unchanged=union(all.shifts, untouchable)
	} else {
		untouchable=NULL
		remain.unchanged=list()
	}
	
	marker=match(s, names(new.vv))
	nn=length(vv)
	K=sum(bb)+1
	N=Ntip(phy)
	
	ncat=sum(bb)
	cur.vv=as.numeric(vv[marker])
	ca.vv=length(which(vv==cur.vv))
	
	if(sum(new.bb)>sum(bb)) {			# add transition: SPLIT
#		print("s")

		n.desc=sum(!all.D%in%remain.unchanged)+1
		n.split=sum(vv==cur.vv)-n.desc
		if(runif(1)<0.5) {
			u=runif(1,min=-0.5*n.desc*cur.vv, max=0.5*n.split*cur.vv)
			nr.desc=cur.vv + u/n.desc
			nr.split=cur.vv - u/n.split
#			(n.desc*nr.desc+n.split*nr.split)/(n.desc+n.split); nr.desc; nr.split
		} else {
			u=runif(1,min=-0.5*n.split*cur.vv, max=0.5*n.desc*cur.vv)
			nr.desc=cur.vv - u/n.desc
			nr.split=cur.vv + u/n.split
#			(n.desc*nr.desc+n.split*nr.split)/(n.desc+n.split); nr.desc; nr.split
		}
		new.vv[vv==cur.vv]=nr.split
		if(length(remain.unchanged)==0) {	# assign new value to all descendants
			new.vv[match(all.D,names(new.vv))] = nr.desc
		} else {							# assign new value to all 'open' descendants 
			new.vv[match(all.D[!(all.D%in%remain.unchanged)],names(new.vv))] = nr.desc
		}
		new.vv[match(s, names(new.vv))]=nr.desc
		lnHastingsRatio = log((K+1)/(2*N-2-K)) ### from Drummond and Suchard 2010: where N is tips, K is number of local parms in tree
		lnPriorRatio = log(ptpois(K+1,lambda,nn)/ptpois(K,lambda,nn))
		
	} else {							# drop transition: MERGE
#		print("m")
		
		anc = get.ancestor.of.node(s, phy)
		if(!is.root(anc, phy)) {			# base new rate on ancestral rate of selected branch
			anc.vv=as.numeric(vv[match(anc,names(vv))])
			na.vv=length(which(vv==anc.vv))
			nr=(anc.vv*na.vv+cur.vv*ca.vv)/(ca.vv+na.vv)
			new.vv[vv==cur.vv | vv==anc.vv]=nr
		} else {							# if ancestor of selected node is root, base new rate on sister node
			sister.tmp=get.desc.of.node(anc,phy)
			sister=sister.tmp[sister.tmp!=s]
			sis.vv=as.numeric(vv[match(sister,names(vv))])
			ns.vv=length(which(vv==sis.vv))
			nr=(sis.vv*ns.vv+cur.vv*ca.vv)/(ca.vv+ns.vv)
			new.vv[vv==cur.vv | vv==sis.vv]=nr			
		}
		lnHastingsRatio = log((2*N-2-K+1)/K) ### from Drummond and Suchard 2010: where N is tips, K is number of local parms in tree
		lnPriorRatio = log(ptpois(K-1,lambda,nn)/ptpois(K,lambda,nn))
		
	}
	
	return(list(new.delta=new.bb, new.values=new.vv, lnHastingsRatio=lnHastingsRatio, lnPriorRatio=lnPriorRatio))
}

switcheroo <- function(cur.values) { 
# mod 11.17.2010 JM Eastman
#	cur.delta=cur.delta.rates; cur.values=cur.rates; theta=FALSE; lambda=log(2); jumpsize=2

	vv=cur.values
	rates.sample=unique(vv)[sample(1:length(unique(vv)), 2)]
	new.vv=vv
	new.vv[which(vv==rates.sample[1])]=rates.sample[2]
	new.vv[which(vv==rates.sample[2])]=rates.sample[1]
	
	return(list(new.values=new.vv))
}

generate.starting.point <- function(data, phy, node.des, theta=FALSE, K=FALSE, jumpsize) { 
# updated JME 11.14.2010
	nn=length(phy$edge.length)
	ntip=Ntip(phy)
	if(!K) nshifts=rtpois(1,log(ntip),nn) else nshifts=K-1
	if(nshifts!=0) bb=sample(c(rep(1,nshifts),rep(0,nn-nshifts)),replace=FALSE) else bb=rep(0,nn)
	names(bb)=phy$edge[,2]
	shifts=as.numeric(names(bb)[bb==1])
	if(theta) {
		min.max=c(adjust.value(min(data), jumpsize), adjust.value(max(data), jumpsize))
		values=runif(sum(bb)+1, min=min.max[min.max==min(min.max)], max=min.max[min.max==max(min.max)])
	} else {
		init.rate=fitContinuous(phy,data)[[1]]$beta
		min.max=c(adjust.rate(init.rate, jumpsize), adjust.rate(init.rate, jumpsize))
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

assign.root.theta <- function(cur.theta, phy) {
	root=phy$edge[1,1]
	desc=phy$edge[which(phy$edge[,1]==root),2]
	root.theta=cur.theta[sample(which(phy$edge[,2]==desc),1)]
	root.theta
}

get.descendants.of.node <- function(node, phy) {
	storage=list()
	if(is.null(node)) node=phy$edge[1,1]
	if(node%in%phy$edge[,1]) {
		desc=c(sapply(node, function(x) {return(phy$edge[which(phy$edge[,1]==x),2])}))
		while(any(desc%in%phy$edge[,1])) {
			desc.true=desc[which(desc%in%phy$edge[,1])]
			desc.desc=lapply(desc.true, function(x) return(get.desc.of.node(x, phy)))
			
			storage=list(storage,union(desc,desc.desc))
			
			desc=unlist(desc.desc)
		}
		storage=sort(unique(unlist(union(desc,storage))))
	}
	return(storage)
}

get.desc.of.node <- function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}

get.ancestor.of.node<-function(node, phy) {
	anc=phy$edge[which(phy$edge[,2]==node),1]
	anc
}

choose.one <- function(cur.delta)
{
	bb=cur.delta
	s=sample(1:length(bb),1)
	bb[s]=1-bb[s]
	bb
}

is.root<-function(node,phy) {
	if(node==phy$edge[1,1]) return(TRUE) else return(FALSE)
}

assess.lnR <- function(lnR) {
	if(is.na(lnR) || abs(lnR)==Inf) {
		error=TRUE
		r=0
	} else {
		if(lnR < -20) {
			r=0 
		} else if(lnR > 0) {
			r=1 
		} else {
			r=exp(lnR)
		}
		error=FALSE
	}
	return(list(r=r, error=error))
}

apply.BMM.to.apetree<-function(apetree, rates) {
	apetree$edge.length=apetree$edge.length*rates
	return(apetree)
}

adjust.value <- function(value, jumpsize) {
# mod 10.20.2010 JM Eastman
	vv=value
	rch <- runif(1, min = -jumpsize, max = jumpsize)
	return(vv[]+rch)
}

adjust.rate <- function(rate, jumpsize) {
# mod 10.20.2010 JM Eastman
	vv=rate
	v=log(vv)
	rch <- runif(1, min = -jumpsize, max = jumpsize)
	v=exp(v[]+rch)
	v
}

prepare.data <- function(phy, data) {
	td <- treedata(phy, data, sort = T)	
	return(list(ape.tre=td$phy, orig.dat=td$data[,1]))
}

summarize.run<-function(cOK, endparms, cur.model) {
#	end = c(nRateProp, nRateSwapProp, nRcatProp, nRootProp, nAlphaProp, nThetaProp, nThetaSwapProp, nTcatProp)

	if(cur.model=="BM") {
		names(endparms)<-names(cOK)<-c("rates", "rate.swap", "rate.catg", "root")
	} else {
		names(endparms)<-names(cOK)<-c("rates", "rate.swap", "rate.catg", "root", "alpha", "theta", "theta.swap", "theta.catg")
	}
	names(endparms)=paste("pr.",names(endparms),sep="")
	names(cOK)=paste("ok.",names(cOK),sep="")
	cat("\n")
	if(cur.model=="BM") {
		print(endparms[paste("pr.", c("rates", "rate.swap", "rate.catg", "root"),sep="")])
		print(cOK[paste("ok.", c("rates", "rate.swap", "rate.catg", "root"),sep="")])
	} else {
		print(endparms[paste("pr.", c("rates", "rate.swap", "rate.catg", "alpha", "theta", "theta.swap", "theta.catg"),sep="")])
		print(cOK[paste("ok.", c("rates", "rate.swap", "rate.catg", "alpha", "theta", "theta.swap", "theta.catg"),sep="")])
	}
}

generate.log<-function(bundled.parms, cur.model, file, init=FALSE) {
#	bundled.parms=list(gen=i, mrate=mean(cur.rates), cats=sum(cur.delta.rates)+1, root=cur.root, alpha=cur.alpha, reg=sum(cur.delta.theta)+1, theta=mean(cur.theta), lnL=curr.lnL)
	res=bundled.parms

	if(init) {
		if(cur.model=="BM") msg=paste("gen", "model", "mean.rate", "rates", "root", "lnL", sep="\t") else msg=paste("gen", "model", "mean.rate", "rates", "alpha", "regimes", "mean.theta", "lnL", sep="\t")
	} else {
		if(cur.model=="BM") {
			msg<-paste(res$gen, sQuote(cur.model), sprintf("%.3f", res$mrate), res$cats, sprintf("%.3f", res$root), sprintf("%.3f", res$lnL),sep="\t")
		} else {
			msg<-paste(res$gen, sQuote(cur.model), sprintf("%.3f", res$mrate), res$cats, sprintf("%.4f", res$alpha), res$reg, sprintf("%.3f", res$theta), sprintf("%.3f", res$lnL),sep="\t")
		}
	}
	write(msg, file=file, append=TRUE)
}

parlogger<-function(model, init=FALSE, node.des, i, curr.lnL, cur.root=NULL, cur.rates=NULL, cur.delta.rates=NULL, cur.alpha=NULL, cur.theta=NULL, cur.delta.theta=NULL, parmBase) {	
	if(init) {
		msg.long=names(node.des[-1])
		if(model=="BM") {
			msg.short=paste("gen", "model", "lnL", sep="\t") 
			parlogs=paste(parmBase, paste(c("summary", "root", "rates", "rate.shifts"), "txt", sep="."), sep="")
			sapply(parlogs[1], function(x) write.table(msg.short, x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			sapply(parlogs[2], function(x) write.table(NULL, x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			sapply(parlogs[3:length(parlogs)], function(x) write.table(paste(msg.long,collapse=" "), x, quote=FALSE, col.names=FALSE, row.names=FALSE))
		} else {
			msg.short=paste("gen", "model", "lnL", sep="\t")
			parlogs=paste(parmBase, paste(c("summary", "root", "alpha", "rates", "rate.shifts", "optima", "optima.shifts"), "txt", sep="."), sep="")
			sapply(parlogs[1], function(x) write.table(msg.short, x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			sapply(parlogs[2:3], function(x) write.table(NULL, x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			sapply(parlogs[4:length(parlogs)], function(x) write.table(paste(msg.long,collapse=" "), x, quote=FALSE, col.names=FALSE, row.names=FALSE))
		}
	} else {
		if(model=="BM") {
			parlogs=paste(parmBase, paste(c("summary", "root", "rates", "rate.shifts"), "txt", sep="."), sep="")
			outlist=list(paste(i, model, curr.lnL, sep="\t"), cur.root, cur.rates, cur.delta.rates) 
			sapply(1:length(outlist), function(x) write.table(paste(outlist[[x]],collapse=" "), parlogs[[x]], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)) 
		} else {
			parlogs=paste(parmBase, paste(c("summary", "root", "alpha", "rates", "rate.shifts", "optima", "optima.shifts"), "txt", sep="."), sep="")
			outlist=list(paste(i, model, curr.lnL, sep="\t"), cur.root, cur.alpha, cur.rates, cur.delta.rates, cur.theta, cur.delta.theta) 
			sapply(1:length(outlist), function(x) write.table(paste(outlist[[x]],collapse=" "), parlogs[[x]], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)) 
		}
	}
}

generate.error.message<-function(i, mod.cur, mod.new, lnR, errorLog, init=FALSE) {
	if(init) {
		initmsg = write(paste("gen", "curr.lnL", "new.lnL", "lnR", sep="\t"), file=errorLog)
	} else {
		write(paste(i, sprintf("%.3f", mod.cur$lnL), sprintf("%.3f", mod.new$lnL), sprintf("%.3f", lnR), sep="\t"), file=errorLog)
	}
}

determine.accepted.proposal<-function(startparms, endparms, cOK) {
	which(startparms!=endparms)->marker
	cOK[marker]=cOK[marker]+1
	cOK
}

rtpois<-function(N, lambda, k) {
	p=ppois(k, lambda, lower.tail=FALSE)
	out=qpois(runif(N, min=0, max=1-p), lambda)
	out
}
				 
ptpois<-function(x, lambda, k) {
	p.k=ppois(k, lambda, lower.tail=FALSE)
	p.x=ppois(x, lambda, lower.tail=FALSE)
	ptp=p.x/(1-p.k)
	return(ptp)
}
				 
## SUMMARIZING MCMC -- PLOTTING FUNCTIONS ##
shifts.plot=function(phy, base.dir, burnin=0, level=0.03, internal.only=FALSE, paint.branches=TRUE, pdf.base=NULL, verbosity.base=NULL, ...) {
	color.length=21
	require(ape)
	oldwd=getwd()
	setwd(base.dir)
	files=dir()
#	if(optima) {
#		shifts.file=files[grep("optima.shifts.txt", files)]
#		estimates.file=files[grep("optima.txt", files)]
#	} else {
		shifts.file=files[grep("rate.shifts.txt", files)]
		estimates.file=files[grep("rates.txt", files)]
#	}
	if(!is.null(pdf.base)) lab=paste(pdf.base, paste("bi", burnin, sep=""), paste("fq", level, sep=""), sep=".")

	cat("READING shifts...\n")
	results=read.table(shifts.file)
	names(results)=as.numeric(results[1,])
	results=results[-c(1:(burnin+1)),]
	if(!all(names(results)%in%phy$edge[,2])) stop("Cannot process results: check that tree matches posterior summaries")
	hits.tmp=apply(results,1,sum)
	hits=length(hits.tmp[hits.tmp>0])
	aa.tmp=apply(results, 2, function(x) sum(x))
	aa.tmp=aa.tmp/hits
	aa=aa.tmp[aa.tmp>=level]
	aa=aa[order(aa, decreasing=TRUE)]
	aa.nodes=aa[which(as.numeric(names(aa))>Ntip(phy))]
	aa.tips=aa[which(as.numeric(names(aa))<=Ntip(phy))]
	aa=aa.nodes
	print(aa.tips)
	nodes=as.numeric(names(aa))
	nodes=nodes[nodes>Ntip(phy)]
	desc=lapply(nodes, function(x) {
				foo=get.descendants.of.node(x,phy)
				if(length(foo)) return(phy$tip.label[foo[foo<=Ntip(phy)]]) else return(NULL)
				}
				)
	shifts=apply(results, 1, sum)
	cat("READING estimates...\n")
	ests=read.table(estimates.file)
	names(ests)=ests[1,]
	ests=ests[-c(1:(burnin+1)),]
	average.ests.tmp=apply(ests,2,mean)
	average.ests=average.ests.tmp-median(average.ests.tmp)
	if(paint.branches) {
		require(colorspace)
#		cce=diverge_hcl(color.length, c = c(100, 0), l = c(50, 90), power = 1.1)
		cce=diverge_hcl(2*color.length, power = 0.5)
		e.seq=seq(min(average.ests),max(average.ests),length=color.length)
		color.gamut=seq(-max(abs(e.seq)), max(abs(e.seq)), length=length(cce))
		colors.branches=sapply(average.ests, function(x) cce[which(min(abs(color.gamut-x))==abs(color.gamut-x))])
		colors.branches=colors.branches[match(as.numeric(names(colors.branches)), phy$edge[,2])]
	} else {
		colors.branches=1
	}
	if(length(nodes)) {
		require(colorspace)
		cols=sapply(nodes, function(x) {
						 a=get.ancestor.of.node(x, phy)
						 comp=ests[,which(names(ests)==a)]
						 this=ests[,which(names(ests)==x)]
						 if(!length(comp)) { # dealing with first descendant of root
						 d=get.desc.of.node(a, phy)
						 d=d[which(d!=x)]
						 
						# find if sister descendant also experienced rate shift 
						 d.shifts=results[, which(names(results)==d)]
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
						 } else {
						 yy=this-comp
						 zz=yy
						 zz[yy>0]=1
						 zz[yy<0]=-1
						 zz[yy==0]=0
						 x.res=mean(zz[zz!=0])
						 }
						 return(x.res)
						 }
						 )
		ccx=diverge_hcl(21, c = c(100, 0), l = c(50, 90), power = 1.1)
		c.seq=round(seq(-1,1,by=0.1),digits=1)
		colors.nodes=ccx[match(round(cols,1),c.seq)]
		names(colors.nodes)=nodes
		colors.nodes=colors.nodes[which(as.numeric(names(colors.nodes))>Ntip(phy))]
	} else {
		colors.nodes=NULL
		cols=NULL
	}
		
# send to pdf 
	if(!is.null(pdf.base)) pdf(file=paste(oldwd, paste(lab,"pdf",sep="."),sep="/"))
	plot(phy, cex=0.1, lwd=0.4, edge.width=0.5, edge.color=colors.branches, ...)
	nn=(Ntip(phy)+1):max(phy$edge[,2])
	ll<-cc<-rr<-rep(0,length(nn))
	ll[nn%in%nodes]=1
	if(length(nodes)) cc[nn%in%nodes]=aa[match(nn[nn%in%nodes],nodes)]/max(aa) 
	if(is.null(colors.nodes)) {
		nodelabels(pch=ifelse(ll==1, 21, NA), bg=ifelse(ll==1, "black", "white"), cex=2*cc)
	} else {
		rr[nn%in%nodes]=colors.nodes[order(names(colors.nodes))]
		nodelabels(pch=ifelse(ll==1, 21, NA), cex=2*cc, col=rr, bg=rr, lwd=0.5)
		if(length(aa.tips) & !internal.only) {
			tt=rep(0,Ntip(phy))
			tt[(1:Ntip(phy))[as.numeric(names(aa.tips))]]=1
			tt.cex=ifelse(tt==1, aa.tips/max(aa.tmp), 0)
			tt.col=sapply(names(aa.tips), function(x) {
						a=get.ancestor.of.node(x, phy)
						comp=ests[,which(names(ests)==a)]
						this=ests[,which(names(ests)==x)]
						if(!length(comp)) { # dealing with first descendant of root
						d=get.desc.of.node(a, phy)
						d=d[which(d!=x)]
						
						# find if sister descendant also experienced rate shift 
						d.shifts=results[, which(names(results)==d)]
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
						} else {
						yy=this-comp
						zz=yy
						zz[yy>0]=1
						zz[yy<0]=-1
						zz[yy==0]=0
						x.res=mean(zz[zz!=0])
						}
						return(x.res)
						}
						)
			tt.ccx=diverge_hcl(21, c = c(100, 0), l = c(50, 90), power = 1.1)
			tt.c.seq=round(seq(-1,1,by=0.1),digits=1)
			tt.colors.nodes=tt.ccx[match(round(tt.col,1),tt.c.seq)]
			names(tt.colors.nodes)=as.numeric(names(aa.tips))
			colors.tips.tmp=tt.colors.nodes[order(as.numeric(names(tt.colors.nodes)))]
			colors.tips=rep(NA, Ntip(phy))
			colors.tips[(1:Ntip(phy))[as.numeric(names(aa.tips))]]=colors.tips.tmp

			tiplabels(pch=ifelse(tt==1, 21, NA), cex=2*tt.cex, col=colors.tips, bg=colors.tips, lwd=0.5)
		}
		legend("bottomleft", title=toupper("direction"), cex=0.5, pt.cex=1, text.col="darkgray", legend = sprintf("%6.1f", rev(c.seq[seq(1,length(c.seq),by=2)])), pch=21, ncol=1, col = "darkgray", pt.bg = rev(ccx[seq(1,length(c.seq),by=2)]), box.lty="blank", border="white")
		if(length(nodes)) legend("bottomright", title=toupper("frequency"), cex=0.5, pt.cex=2*(seq(min(cc), max(cc), length=6)), text.col="darkgray", legend = sprintf("%10.3f", seq(round(min(aa),2), round(max(aa),2), length=6)), pch=21, ncol=1, col = "darkgray", pt.bg = "white", box.lty="blank", border="white")
	}
	if(!is.null(pdf.base)) dev.off()
		
	if(length(nodes)) names(cols)=names(aa)
	setwd(oldwd)
	names(desc)=paste("node",nodes,sep=".")
	if(!is.null(verbosity.base)) {
		summary.table=c(hits/nrow(results), mean(shifts), median(shifts))
		names(summary.table)=c("shift.prob", "mean.shift","median.shifts")
		write.table(summary.table, file=paste(verbosity.base, "run.summary.txt", sep="."), quote=FALSE, col.names=FALSE)
		if(length(nodes)) {
			for(nn in 1:length(nodes)) {
				write.table(c(nodes[nn], desc[[nn]]), file=paste(verbosity.base, "shift.descendants.txt", sep="."), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
			}
		}
		allres=merge(as.data.frame(aa),as.data.frame(cols),by=0)
		allres=allres[order(allres$aa, decreasing=TRUE),]
		names(allres)=c("node", "conditional.shift.probability", "relative.shift.direction")
		write.table(allres, file=paste(verbosity.base, "estimates.summary.txt", sep="."), quote=FALSE, row.names=FALSE)
	}
	return(list(desc=desc, shift.prob=hits/nrow(results), mean.shifts=mean(shifts), median.shifts=median(shifts), shift.prop=aa, shift.direction=cols, mean.rates=average.ests.tmp, mean.shifts=aa.tmp/nrow(results)))
}


plot.quant<-function(phy, data, data.names, shrinker=1, relative=TRUE, cex=0.5, adj=c(15,5), lwd=0.5, font=1, ...){
	std=function(x,max,min){return((x-min)/(max-min))}
	if(relative) {
		data=sapply(data, function(x) std(x, max(data), min(data)))
	}
	
	f=match(phy$tip.label, data.names)
	if(is.na(sum(f))) {
		stop("tips cannot be matched to tree")
	}
	else {
		data=data[f]
		phy$trait=data
		plot.phylo(phy, show.tip.label=T, label.offset=adj[1], cex=cex, font=font, ...)
		tiplabels(text=NULL, pch=21, adj=adj[2], cex=c(phy$trait/shrinker), bg="white",lwd=lwd)
		} 
}

branchcol.plot=function(phy, cur.rates, color.length=4, ...) {
	color.length=color.length
	require(ape)
	require(colorspace)
	ests=cur.rates
	names(ests)=phy$edge[,2]
	cce=diverge_hcl(color.length, power = 0.5)
	e.seq=seq(min(ests),max(ests),length=color.length)
	colors.branches=sapply(ests, function(x) cce[which(min(abs(e.seq-x))==abs(e.seq-x))])
	names(colors.branches)=names(ests)
	colors.branches=colors.branches[match(as.numeric(names(colors.branches)), phy$edge[,2])]	
	plot(phy, cex=0.1, lwd=0.4, edge.width=0.5, edge.color=colors.branches, ...)
}


prior.posterior.plot=function(data, rpois.lambda) {
# data should be a vector of number of shifts, for instance, for all interations across the MCMC sampling (post-burnin)
	require(ggplot2)
	dataset <- data.frame(variable = gl(2, length(data), labels = c("posterior", "prior")), value = c(data, rpois(length(data),rpois.lambda))) 
	ggplot(data = dataset, aes(x = value, fill = variable)) + geom_histogram(position = "dodge") + scale_fill_manual(values = c("gray", "black"))	
}



