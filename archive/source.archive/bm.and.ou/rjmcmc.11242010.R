## Written by LJ HARMON (2008), AL HIPP (2008), and JM EASTMAN (2010)
# MODELING TRAIT EVOLUTION by BROWNIAN MOTION

rjMCMC.bm<-function (phy, data, ngen=1000, sampleFreq=100, probMergeSplit=0.05, probRoot=0.01, lambdaK=log(2), constrainK=FALSE, jumpsize=2, fileBase="result", simplestart=FALSE) 
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
			init.rate	<- generate.starting.point(orig.dat, cur.tre, node.des, theta=FALSE, K=constrainK, jumpsize=jumpsize)
		} else {
			init.rate	<- list(values=rep(fitContinuous(cur.tre,orig.dat)$Trait1$beta,length(cur.tre$edge.length)),delta=rep(0,length(cur.tre$edge.length)))
		}
	} else {
		init.rate	<- generate.starting.point(orig.dat, cur.tre, node.des, theta=FALSE, K=constrainK, jumpsize=jumpsize )
	}
	
    cur.rates		<- init.rate$values			
	cur.delta.rates	<- init.rate$delta
	cur.root		<- adjust.value(mean(orig.dat), jumpsize)
	cur.tre			<- apply.BMM.to.apetree(cur.tre, cur.rates)
	
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
			cur.tre <- new.tre
			mod.cur <- mod.new
			curr.lnL <- mod.new$lnL
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

rjMCMC.ou<-function (phy, data, ngen=1000, sampleFreq=100, probMergeSplit=0.05, probRoot=0.01, probRate=0.1, probAlpha=0.01, upperAlpha=100, lambdaK=log(2), constrainK=FALSE, jumpsize=2, fileBase="result", simplestart=FALSE) 
{
#	load("rjMCMC/versions/fakedata.rda"); phy=sims$phy; data=sims[[2]]; ngen=1000; sampleFreq=100; probMergeSplit=0.05; probRoot=0.01; probRate=0.1; probAlpha=0.01; lambdaK=log(2); constrainK=FALSE; jumpsize=2; upperAlpha=100; fileBase="result"; simplestart=FALSE 
	
	model="OU"
	heat=1
	require(ouch)
	require(geiger)
	abs.low.alpha=1e-13
	if(length(upperAlpha)==2 & upperAlpha[1]>abs.low.alpha) {
		alpha.bounds=upperAlpha
	} else if(!is.null(upperAlpha) & is.numeric(upperAlpha)){
		alpha.bounds=c(abs.low.alpha, upperAlpha)
	} else {
		alpha.bounds=c(abs.low.alpha, 100)
	}
	
### prepare data for rjMCMC
	cur.model		<- model								
	
	dataList		<- prepare.ouchdata(phy, data)			
	ape.tre			<- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat	
	ouch.tre		<- dataList$ouch.tre
	ouch.dat		<- dataList$ouch.dat
	ape.to.ouch.order <- match(as.numeric(names(dataList$ape.ordered[-1])), as.numeric(ape.tre$edge[,2]))
	nn				<- length(ape.tre$edge.length)	
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))
	subtrees		<- subtrees(phy)
	subtrees.marker <- sapply(subtrees, function(x) return(min(x$node.label)))
	
# initialize parameters
	if(is.numeric(constrainK) & (constrainK > length(ape.tre$edge) | constrainK < 1)) stop("Constraint on rate shifts is nonsensical. Ensure that constrainK is at least 1 and less than the number of available nodes in the tree.")
	
	if(simplestart | is.numeric(constrainK) ) {
		if(is.numeric(constrainK)) {
			init.theta	<- generate.starting.point(orig.dat, ape.tre, node.des, theta=TRUE, K=constrainK, jumpsize=jumpsize)
		} else {
			init.theta	<- list(values=rep(mean(orig.dat),length(ape.tre$edge.length)),delta=rep(0,length(ape.tre$edge.length)))
		}
	} else {
		init.theta	<- generate.starting.point(orig.dat, ape.tre, node.des, theta=TRUE, K=constrainK, jumpsize=jumpsize )
	}
	
	cur.rate		<- fitContinuous(ape.tre,orig.dat)$Trait1$beta
	cur.alpha		<- min(runif(1000, alpha.bounds[1], alpha.bounds[2]))
    cur.theta		<- init.theta$values			
	cur.delta.theta	<- init.theta$delta
	cur.root.theta	<- assign.root.theta(cur.theta,ape.tre)
	cur.regimes		<- c(cur.root.theta,cur.theta) 
	
	mod.cur = new.ou.lik.fn(cur.rate, cur.alpha, cur.regimes, cur.theta, ouch.dat, ouch.tre, ape.to.ouch.order)
	
# proposal counts
	nRateProp =		0	
	nRootProp =		0
	nAlphaProp =	0										
	nTcatProp =		0										
	nThetaProp =	0
	nRateOK =		0											
	nRootOK =		0
	nAlphaOK =		0										
	nTcatOK =		0
	nThetaOK =		0
	
	cOK=c(												
		  nRateOK,
		  nRootOK,
		  nAlphaOK,
		  nTcatOK,
		  nThetaOK
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
	
	probSwap=max(c(probRoot, probRate, probAlpha))
	minorprobs=sum(c(probRoot, probRate, probAlpha, probSwap))
### Begin rjMCMC
    for (i in 1:ngen) {
		
		# loop to ensure that each state produces a sensible proposal 
		while(1) {
			lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
			startparms = c(nRateProp, nRootProp, nAlphaProp, nTcatProp, nThetaProp)
			
			## OU IMPLEMENTATION ##
			if (runif(1) < (2 * probMergeSplit) & !constrainK) {									# adjust theta categories
				nr=split.or.merge(cur.delta.theta, cur.theta, ape.tre, node.des, lambdaK)
				new.alpha=cur.alpha
				new.theta=nr$new.values ##
				new.delta.theta=nr$new.delta ##
				new.root.theta<-nrt<-assign.root.theta(new.theta, ape.tre) ## 
				new.regimes=c(nrt,new.theta) ##
				new.rate=cur.rate 
				nTcatProp=nTcatProp+1
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
				
			} else { # neither split nor merge
				while(1) {
					if(runif(1)<(probAlpha/minorprobs)) {											# adjust alpha
						new.alpha=adjust.rate(cur.alpha, jumpsize) ##
						new.root.theta=cur.root.theta
						new.theta=cur.theta
						new.delta.theta=cur.delta.theta
						new.regimes=cur.regimes 
						new.rate=cur.rate
						nAlphaProp=nAlphaProp+1
						break()
					} else if(runif(1)<(probRoot/minorprobs)) {										# adjust root
						new.alpha=cur.alpha									
						new.root.theta=assign.root.theta(cur.theta, ape.tre) ##
						new.theta=cur.theta
						new.delta.theta=cur.delta.theta
						new.regimes=cur.regimes 
						new.rate=cur.rate 
						nRootProp=nRootProp+1					
						break()
					} else if(runif(1)<(probRate/minorprobs)) {										# adjust rate
						new.alpha=cur.alpha									
						new.root.theta=cur.root.theta
						new.theta=cur.theta
						new.delta.theta=cur.delta.theta
						new.regimes=cur.regimes 
						new.rate=adjust.rate(cur.rate, jumpsize) ##
						nRateProp=nRateProp+1
						break()
					} else if(runif(1)<0.5*(probSwap/minorprobs) & length(unique(cur.theta))>1) {	# swap thetas
						nr=switcheroo(cur.theta)
						new.alpha=cur.alpha 
						new.theta=nr$new.values ##
						new.delta.theta=cur.delta.theta 
						new.root.theta<-nrt<-assign.root.theta(new.theta, ape.tre) ## 
						new.regimes=c(nrt,new.theta) ##
						new.rate=cur.rate
						nThetaProp=nThetaProp+1					
						break()
					} else {																		# adjust thetas
						new.alpha=cur.alpha 
						new.theta=adjust.value(cur.theta, jumpsize, theta=TRUE) ##
						new.delta.theta=cur.delta.theta 
						new.root.theta<-nrt<-assign.root.theta(new.theta, ape.tre) ## 
						new.regimes=c(nrt,new.theta) ##
						new.rate=cur.rate
						nThetaProp=nThetaProp+1					
						break()
					}
					
				}
			}
			if(within.range(new.alpha, alpha.bounds[1], alpha.bounds[2])) {
				mod.new=try(new.ou.lik.fn(new.rate, new.alpha, new.regimes, new.theta, ouch.dat, ouch.tre, ape.to.ouch.order),silent=TRUE)
				if(!inherits(mod.new, "try-error")) break()
			}

		}
		lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
	
		
# compare likelihoods
		endparms = c(nRateProp, nRootProp, nAlphaProp, nTcatProp, nThetaProp)

		r=assess.lnR((heat * lnLikelihoodRatio + heat * lnPriorRatio + lnHastingsRatio)->lnR)
		
		if (runif(1) <= r$r) {			# adopt proposal
			decision="adopt"
			cur.alpha <- new.alpha
			cur.theta <- new.theta
			cur.delta.theta <- new.delta.theta
			cur.root.theta <- new.root.theta
			cur.regimes <- new.regimes
			cur.rate <- new.rate
			cOK <- determine.accepted.proposal(startparms, endparms, cOK)
			mod.cur <- mod.new
			curr.lnL <- mod.new$lnL
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
			bundled.parms=list(gen=i, rate=cur.rate, regimes=sum(cur.delta.theta)+1, root=cur.root.theta, mean.theta=mean(cur.theta), alpha=cur.alpha, lnL=curr.lnL)
			generate.log(bundled.parms, cur.model, file=runLog)			
			parlogger(model=model, init=FALSE, node.des=node.des, i=i, curr.lnL=curr.lnL, cur.root=cur.root.theta, cur.rates=cur.rate, cur.alpha=cur.alpha, cur.theta=cur.theta, cur.delta.theta=cur.delta.theta, parmBase=parmBase)
		}
	}
	
# End rjMCMC
	close(errorLog)
	close(runLog)
	summarize.run(cOK, endparms, cur.model)	
}


## AUXILIARY FUNCTIONS

new.ou.lik.fn <- function(rate, alpha, regimes, theta, ouch.dat, ouch.tre, ape.to.ouch.order) { # from OUCH 2.7-1:::hansen()
# mod 11.17.2010 JM Eastman; checked -- returns very similar if tested against fitted values from hansen()
#	rate=new.rate; alpha=new.alpha; regimes=new.regimes; theta=new.theta
	ouch.dat$regimes=as.factor(regimes[c(1,ape.to.ouch.order+1)])
#	theta=sort(unique(theta))
	alpha.ouch=as.matrix(alpha)
	n <- length(which(!is.na(ouch.dat['trait'])))	
	beta=ouch:::regime.spec(ouch.tre, ouch.dat['regimes']) 
	theta=as.numeric(levels(ouch:::sets.of.regimes(ouch.tre,ouch.dat['regimes'])$regimes))
	dat=ouch.dat$trait[n:nrow(ouch.dat)]
	names(dat)=ouch.dat$labels[n:nrow(ouch.dat)]
	ev <- eigen(alpha.ouch,symmetric=TRUE)
	w <- .Call(ouch:::ouch_weights,object=ouch.tre,lambda=ev$values,S=ev$vectors,beta=beta) 
	v <- .Call(ouch:::ouch_covar,object=ouch.tre,lambda=ev$values,S=ev$vectors,sigma.sq=rate)
	y <- matrix(theta,ncol=1)
	e <- w%*%y-dat
	dim(y) <- dim(y)[1]
	dim(e) <- n
	q <- e%*%solve(v,e)
	det.v <- determinant(v,logarithm=TRUE)
	dev <- n*log(2*pi)+as.numeric(det.v$modulus)+q[1,1]
	list(
		 theta=theta,
		 weight=w,
		 vcov=v,
		 resids=e,
		 lnL=-0.5*dev
		 )
}

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
		decision="split"
		n.desc=sum(!all.D%in%remain.unchanged)+1
		n.split=sum(vv==cur.vv)-n.desc
		if(runif(1)<0.5) {
			uu=c(-0.5*n.desc*cur.vv, 0.5*n.split*cur.vv)
			u=runif(1,min=min(uu), max=max(uu))
			nr.desc=cur.vv + u/n.desc
			nr.split=cur.vv - u/n.split
#			(n.desc*nr.desc+n.split*nr.split)/(n.desc+n.split); nr.desc; nr.split
		} else {
			uu=c(-0.5*n.split*cur.vv, 0.5*n.desc*cur.vv)
			u=runif(1,min=min(uu), max=max(uu))
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
		decision="merge"
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
	
	return(list(new.delta=new.bb, new.values=new.vv, lnHastingsRatio=lnHastingsRatio, lnPriorRatio=lnPriorRatio, decision=decision))
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
			   
within.range=function(x, min, max) {			 
	a=sign(x-min)
	b=sign(x-max)
	if(abs(a+b)==2) return(FALSE) else return(TRUE)
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

ape.to.ouch.tre<-function (tree, scale = FALSE, branch.lengths = tree$edge.length) 
{
    if (!inherits(tree, "phylo")) 
	stop(sQuote("tree"), " must be of class ", sQuote("phylo"))
    nnodes <- nrow(tree$edge) + 1
    n.term <- length(tree$tip.label)
    n.int <- nnodes - n.term
    tmp <- matrix(NA, nnodes, 2)
    tmp[-1, 1:2] <- tree$edge
    tmp[1, 2] <- nnodes + 1
    bl <- c(0, branch.lengths)
	bl.ape=bl
	names(bl.ape)=c("root", tree$edge[,2])
    reord <- order(-tmp[, 2])
    bl <- bl[reord]
	bl.ouch <- bl.ape[reord]
    tmp <- tmp[reord, ]
    node <- seq(nnodes)
    ancestor <- rep(NA, nnodes)
    for (n in 2:nnodes) {
        anc <- which(tmp[, 2] == tmp[n, 1])
        if (length(anc) > 1) 
		stop("invalid tree")
        if (length(anc) > 0) {
            ancestor[n] <- node[anc]
        }
        else {
            ancestor[n] <- node[1]
        }
    }
    if (is.null(tree$node.label)) 
	tree$node.label <- rep("", n.int)
    species <- rev(c(tree$tip.label, tree$node.label[-1], tree$node.label[1]))
    times <- rep(NA, nnodes)
    for (n in 1:nnodes) times[n] <- ouch:::branch.height(node, ancestor, bl, n)
    if (is.logical(scale)) {
        if (is.na(scale)) 
		stop("if ", sQuote("scale"), " is logical, it must be either true or false")
        if (scale) 
		times <- times/max(times)
    }
    else if (is.numeric(scale)) {
        times <- times/abs(scale)
    }
    else {
        stop(sQuote("scale"), " must be either logical or numeric")
    }
    return(list(ouch.tre=ouchtree(nodes = node, ancestors = ancestor, times = times, labels = species), ape.ordered=bl.ouch))
}

adjust.value <- function(value, jumpsize, theta=FALSE) {
# mod 10.20.2010 JM Eastman
	vv=value
	if(runif(1)<0.5 | length(unique(vv))==1) {
		rch <- runif(1, min = -jumpsize, max = jumpsize)
		return(vv[]+rch)
	} else {
		return((vv-mean(vv))*runif(1,min=-jumpsize, max=jumpsize)+mean(vv))
	}
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

prepare.ouchdata <- function(phy, data) {
	td <- treedata(phy, data, sort = T)	
	
	ot <- ape.to.ouch.tre(td$phy->ape.tre, scale=FALSE)		
	ouch.tre=ot$ouch.tre
	ouch.dat=data.frame(as.numeric(ouch.tre@nodes), as.character(ouch.tre@nodelabels))
	names(ouch.dat)=c("nodes","labels")
	ouch.dat[ouch.dat[,2]=="",2]=NA
	ouch.dat$trait=NA
	mM=match(ouch.dat$labels,names(data))
	for(i in 1:nrow(ouch.dat)){
		if(!is.na(mM[i])){
			ouch.dat$trait[i]=data[mM[i]]	
		}
	}
	return(list(ouch.tre=ouch.tre, ouch.dat=ouch.dat, ape.tre=td$phy, orig.dat=td$data[,1], ape.ordered=ot$ape.ordered))
}

summarize.run<-function(cOK, endparms, cur.model) {
#	endparmsOU = c(nRateProp, nRootProp, nAlphaProp, nTcatProp, nThetaProp)

	if(cur.model=="BM") {
		names(endparms)<-names(cOK)<-c("rates", "rate.swap", "rate.catg", "root")
	} else {
		names(endparms)<-names(cOK)<-c("rate", "root", "alpha", "regimes", "theta")
	}
	names(endparms)=paste("pr.",names(endparms),sep="")
	names(cOK)=paste("ok.",names(cOK),sep="")
	cat("\n")
	if(cur.model=="BM") {
		print(endparms[paste("pr.", c("rates", "rate.swap", "rate.catg", "root"),sep="")])
		print(cOK[paste("ok.", c("rates", "rate.swap", "rate.catg", "root"),sep="")])
	} else {
		print(endparms[paste("pr.", c("rate", "root", "alpha", "regimes", "theta"),sep="")])
		print(cOK[paste("ok.", c("rate", "root", "alpha", "regimes", "theta"),sep="")])
	}
}

generate.log<-function(bundled.parms, cur.model, file, init=FALSE) {
#	bundled.parms=list(gen=i, mrate=mean(cur.rates), cats=sum(cur.delta.rates)+1, root=cur.root, alpha=cur.alpha, reg=sum(cur.delta.theta)+1, theta=mean(cur.theta), lnL=curr.lnL)
	res=bundled.parms

	if(init) {
		if(cur.model=="BM") msg=paste("gen", "model", "mean.rate", "rates", "root", "lnL", sep="\t") else msg=paste("gen", "model", "rate", "root", "alpha", "regimes", "mean.theta", "lnL", sep="\t")
	} else {
		if(cur.model=="BM") {
			msg<-paste(res$gen, sQuote(cur.model), sprintf("%.3f", res$mrate), res$cats, sprintf("%.3f", res$root), sprintf("%.3f", res$lnL),sep="\t")
		} else {
			msg<-paste(res$gen, sQuote(cur.model), sprintf("%.3f", res$rate), sprintf("%.2f", res$root), sprintf("%.4f", res$alpha), res$regimes, sprintf("%.3f", res$mean.theta), sprintf("%.3f", res$lnL),sep="\t")
		}
	}
	write(msg, file=file, append=TRUE)
}

parlogger<-function(model, init=FALSE, node.des, i, curr.lnL, cur.root=NULL, cur.rates=NULL, cur.delta.rates=NULL, cur.alpha=NULL, cur.theta=NULL, cur.delta.theta=NULL, parmBase) {	
	if(init) {
		msg.long=names(node.des[-1])
		if(model=="BM") {
			msg.short=paste("gen", "root", "lnL", sep="\t") 
			parlogs=paste(parmBase, paste(c("summary", "rates", "rate.shifts"), "txt", sep="."), sep="")
			sapply(parlogs[1], function(x) write.table(msg.short, x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			sapply(parlogs[2:length(parlogs)], function(x) write.table(paste(msg.long,collapse=" "), x, quote=FALSE, col.names=FALSE, row.names=FALSE))
		} else {
			msg.short=paste("gen", "root", "rate", "alpha", "lnL", sep="\t")
			parlogs=paste(parmBase, paste(c("summary", "optima", "optima.shifts"), "txt", sep="."), sep="")
			sapply(parlogs[1], function(x) write.table(msg.short, x, quote=FALSE, col.names=FALSE, row.names=FALSE))
			sapply(parlogs[2:length(parlogs)], function(x) write.table(paste(msg.long,collapse=" "), x, quote=FALSE, col.names=FALSE, row.names=FALSE))
		}
	} else {
		if(model=="BM") {
			parlogs=paste(parmBase, paste(c("summary", "rates", "rate.shifts"), "txt", sep="."), sep="")
			outlist=list(paste(i, cur.root, curr.lnL, sep="\t"), cur.rates, cur.delta.rates) 
			sapply(1:length(outlist), function(x) write.table(paste(outlist[[x]],collapse=" "), parlogs[[x]], quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)) 
		} else {
			parlogs=paste(parmBase, paste(c("summary", "optima", "optima.shifts"), "txt", sep="."), sep="")
			outlist=list(paste(i, cur.root, cur.rates, cur.alpha, curr.lnL, sep="\t"), cur.theta, cur.delta.theta) 
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
shifts.plot=function(phy, base.dir, burnin=0, level=0.03, internal.only=FALSE, paint.branches=TRUE, pdf=TRUE, verbose=TRUE, ...) {
#	phy; base.dir; burnin=0; level=0.03; internal.only=FALSE; paint.branches=TRUE; pdf=TRUE; verbose=TRUE
	nullphy=phy
	nullphy$tip.label[]=""
	color.length=21
	require(ape)
	oldwd=getwd()
	setwd(base.dir)
	files=dir()
	if(any(grep("optima", files))) optima=TRUE else optima=FALSE
	if(optima) {
		shifts.file=files[grep("optima.shifts.txt", files)]
		estimates.file=files[grep("optima.txt", files)]
	} else {
		shifts.file=files[grep("rate.shifts.txt", files)]
		estimates.file=files[grep("rates.txt", files)]
	}
	if(pdf | verbose) {
		lab=paste(strsplit(base.dir, ".", fixed=TRUE)[[1]][1:2],collapse=".")
		out.base=paste(lab, paste("bi", burnin, sep=""), paste("fq", level, sep=""), sep=".")
		outdir=paste("./results")
		if(!file.exists(outdir)) dir.create(outdir)
	}
	cat("READING shifts...\n")
	results=read.table(shifts.file)
	names(results)=as.numeric(results[1,])
	burnin=ceiling(burnin*nrow(results))
	results=results[-c(1:(burnin+1)),]
	if(!all(names(results)%in%phy$edge[,2])) stop("Cannot process results: check that tree matches posterior summaries")
	hits.tmp<-shifts<-apply(results,1,sum)
	hits=length(hits.tmp[hits.tmp>0])
	aa.tmp=apply(results, 2, function(x) sum(x))
	aa.tmp=aa.tmp/hits
	aa=aa.tmp[aa.tmp>=level]
	aa=aa[order(aa, decreasing=TRUE)]
	aa.nodes=aa[which(as.numeric(names(aa))>Ntip(phy))]
	aa.tips=aa[which(as.numeric(names(aa))<=Ntip(phy))]
	aa=aa.nodes
	nodes=as.numeric(names(aa))
	nodes=nodes[nodes>Ntip(phy)]
	desc=lapply(nodes, function(x) {
				foo=get.descendants.of.node(x,phy)
				if(length(foo)) return(phy$tip.label[foo[foo<=Ntip(phy)]]) else return(NULL)
				}
				)
	cat("READING estimates...\n")
	ests=read.table(estimates.file)
	names(ests)=ests[1,]
	ests=ests[-c(1:(burnin+1)),]
	average.ests.tmp=apply(ests,2,mean)
	average.ests=average.ests.tmp-median(average.ests.tmp)
	if(paint.branches) {
		require(colorspace)
		cce=diverge_hcl(2*color.length, power = 0.5)
		e.seq=seq(min(average.ests),max(average.ests),length=color.length)
		color.gamut=seq(-max(abs(e.seq)), max(abs(e.seq)), length=length(cce))
		colors.branches=sapply(average.ests, function(x) {y=cce[which(min(abs(color.gamut-x))==abs(color.gamut-x))]; return(y[sample(1:length(y),1)])})
		colors.branches=colors.branches[match(as.numeric(names(colors.branches)), phy$edge[,2])]
	} else {
		colors.branches=1
	}
	ccx=diverge_hcl(color.length, c = c(100, 0), l = c(50, 90), power = 1.1)
	c.seq=round(seq(-1,1,by=0.1),digits=1)
	
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
					return(x.res)
					}
					)
		colors.nodes=ccx[match(round(cols,1),c.seq)]
		names(colors.nodes)=nodes
		colors.nodes=colors.nodes[which(as.numeric(names(colors.nodes))>Ntip(phy))]
	} else {
		colors.nodes=NULL
		cols=NULL
	}
	
	if(length(aa.tips) & !internal.only) {
		tt<-tt.cex<-rep(0,Ntip(phy))
		tt[as.numeric(names(aa.tips))]=1
#		tt.cex[as.numeric(names(aa.tips))]=aa.tips/max(aa.tmp)
		tt.cex[as.numeric(names(aa.tips))]=aa.tips
		
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
		colors.tips[as.numeric(names(colors.tips.tmp))]=colors.tips.tmp
		tip.names=as.list(phy$tip.label[as.numeric(names(aa.tips))])
		names(tip.names)=paste("node",names(aa.tips),sep=".")
	}	else {
		tt.cex=NULL
		tt.col=NULL
		tip.names=NULL
		colors.tips=NULL
		tt=rep(0,Ntip(phy))
	}
	
## PLOTTING OF TREE
	CEX=max(c(0.05, (2+-0.36*log(Ntip(phy))) ))	
	if(pdf) pdf(file=paste(outdir, paste(out.base,"pdf",sep="."),sep="/"))
	plot(phy, cex=CEX, lwd=0.4, edge.width=0.5, edge.color=colors.branches, show.tip.label=TRUE, label.offset=strwidth(par("pch"),cex=1.25), ...)
	nn=(Ntip(phy)+1):max(phy$edge[,2])
	ll<-cc<-rr<-rep(0,length(nn))
	ll[nn%in%nodes]=1
#	if(length(nodes)) cc[nn%in%nodes]=aa[match(nn[nn%in%nodes],nodes)]/max(c(aa,tt.cex)) 
	if(length(nodes)) cc[nn%in%nodes]=aa[match(nn[nn%in%nodes],nodes)]
	rr[nn%in%nodes]=colors.nodes[order(names(colors.nodes))]
	nodelabels(pch=ifelse(ll==1, 21, NA), cex=2*cc, col=rr, bg=rr, lwd=0.5)
	
	tiplabels(pch=ifelse(tt==1, 21, NA), cex=2*tt.cex, col=colors.tips, bg=colors.tips, lwd=0.5)
	legend("topright", title=toupper("direction"), cex=0.5, pt.cex=1, text.col="darkgray", legend = sprintf("%6.1f", rev(c.seq[seq(1,color.length,by=2)])), pch=21, ncol=1, col = "darkgray", pt.bg = rev(ccx[seq(1,color.length,by=2)]), box.lty="blank", border="white")
	if(length(nodes) | (length(aa.tips) & !internal.only)) {
		legend("bottomright", title=toupper("frequency"), cex=0.5, pt.cex=2*(seq(0, 1, length=9)), text.col="darkgray", legend = sprintf("%10.3f", seq(0, 1, length=9)), pch=21, ncol=1, col = "darkgray", pt.bg = "white", box.lty="blank", border="white")
	}
	if(pdf) dev.off()
	
	if(length(nodes)) {
		names(cols)=names(aa)
		names(desc)=paste("node",nodes,sep=".")
	}
##
	
	if(verbose) {
		summary.table=c(hits/nrow(results), mean(shifts), median(shifts))
		names(summary.table)=c("shift.prob", "mean.shift","median.shifts")
		write.table(summary.table, file=paste(outdir, paste(out.base, "run.summary.txt", sep="."), sep="/"), quote=FALSE, col.names=FALSE)
		if(length(nodes)) {
			for(nn in 1:length(nodes)) {
				write.table(c(nodes[nn], desc[[nn]]), file=paste(outdir, paste(out.base, "shift.descendants.txt", sep="."), sep="/"), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
			}
		}
		if(length(aa.tips)) {
			for(nn in 1:length(aa.tips)) {
				write.table(c(as.numeric(names(aa.tips)[nn]), tip.names[[nn]]), file=paste(outdir, paste(out.base, "shift.descendants.txt", sep="."), sep="/"), append=TRUE, quote=FALSE, row.names=FALSE, col.names=FALSE)
			}
		}
		allres=data.frame(matrix(NA, nrow=nrow(phy$edge), ncol=3))
		if(length(aa)) {
			x=as.numeric(names(aa))
			nodres=merge(as.data.frame(aa),as.data.frame(cols),by=0)
			allres[nodres[,1],]=nodres
		}
		if(length(aa.tips)) {
			x=as.numeric(names(aa.tips))
			tipres=merge(as.data.frame(aa.tips), as.data.frame(tt.col),by=0)
			allres[tipres[,1],]=tipres
		}
		allres=allres[!apply(allres, 1, function(x)any(is.na(x))),]
		allres=allres[order(allres[,2], decreasing=TRUE),]
		if(nrow(allres)>0) {
			allres=data.frame(allres)
			names(allres)=c("node", "conditional.shift.probability", "relative.shift.direction")
			write.table(allres, file=paste(outdir, paste(out.base, "estimates.summary.txt", sep="."), sep="/"), quote=FALSE, row.names=FALSE)
		}
	}
	setwd(oldwd)
	
	if(!optima) {
		return(list(desc=desc, tips=tip.names, shift.prob=hits/nrow(results), mean.shifts=mean(shifts), median.shifts=median(shifts), nodeshift.prop=aa, nodeshift.direction=cols, tipshift.prop=aa.tips, tipshift.direction=tt.col, mean.rates=average.ests.tmp, mean.shifts=aa.tmp/nrow(results)))
	} else {
		return(list(desc=desc, tips=tip.names, shift.prob=hits/nrow(results), mean.shifts=mean(shifts), median.shifts=median(shifts), nodeshift.prop=aa, nodeshift.direction=cols, tipshift.prop=aa.tips, tipshift.direction=tt.col, mean.optima=average.ests.tmp, mean.shifts=aa.tmp/nrow(results)))
	}
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

prior.posterior.plot.alt=function(outlogs, lambda, ...){
	foo=outlogs
	d=lapply(dir(foo),function(x) read.table(paste(foo,x,sep="/"),header=T))
	names=dir(foo)
	lapply(names,function(x)paste(strsplit(x,".",fixed=T)[[1]][1:2],collapse="."))->nn
	
	
	n=length(nn)
	ss=(1:10)^2
	sn=ss-n
	s=ss[which(sn==min(sn[sn>=0]))]
	x=sqrt(s)
	layout(matrix(1:s, x, x))
		
	for(i in 1:length(d)){
		if(any(names(d[[i]])=="regimes")) {
			plot(density(d[[i]]$regimes), bty="n",xlim=c(0,10), main=nn[[i]],xlab="regimes") 
		} else {
			plot(density(d[[i]]$rates),bty="n",xlim=c(0,10), main=nn[[i]],xlab="rates")
		}
		lines(density(rpois(nrow(d[[i]]),lambda)), ...)
	}	
}

