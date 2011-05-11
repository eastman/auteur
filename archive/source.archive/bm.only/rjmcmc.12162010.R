## Written by LJ HARMON (2008), AL HIPP (2008), and JM EASTMAN (2010)
# MODELING TRAIT EVOLUTION by BROWNIAN MOTION

rjMCMC.bm<-function (phy, data, ngen=1000, sample.freq=100, prob.mergesplit=0.05, prob.root=0.01, lambdaK=log(2), constrainK=FALSE, jumpsize=NULL, simplestart=FALSE, internal.only=FALSE, summary=TRUE, fileBase="result") 
{ 
colorplots=FALSE	
# revisions to code:
	
# 12.15.2010: allowing chain to update jumpsize (with calibration before run as well as during first 10th of chain)  JME
# 12.15.2010: revised proposal mechanism to ensure that proposals are conducted at rates specified by user JME
# 12.15.2010: dropped switcheroo proposal JME
# 12.16.2010: switched split.or.merge again to log-space
	
# ngen=1000; sample.freq=100; prob.mergesplit=0.05; prob.root=0.01; lambdaK=log(2); constrainK=FALSE; jumpsize=NULL; fileBase="result"; simplestart=FALSE; summary=TRUE 

	model="BM"
	heat=1
	require(geiger)
	
	if(is.null(jumpsize)) {
		cat("CALIBRATING jumpsize...\n")
		adjustable.jumpsize=TRUE
		jumpsize=calibrate.jumpsize(phy, data, nsteps=ngen/1000)
	} else {
		adjustable.jumpsize=FALSE
	}
	

### prepare data for rjMCMC
	cur.model		<- model
	dataList		<- prepare.data(phy, data)			
	ape.tre			<- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat						
	nn				<- length(ape.tre$edge.length)	
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))
	subtrees		<- subtrees(phy)
	subtrees.marker <- sapply(subtrees, function(x) return(min(x$node.label)))

# initialize parameters
	if(is.numeric(constrainK) & (constrainK > length(ape.tre$edge) | constrainK < 1)) stop("Constraint on rate shifts is nonsensical. Ensure that constrainK is at least 1 and less than the number of available nodes in the tree.")

	if(simplestart | is.numeric(constrainK) | internal.only) {
		if(is.numeric(constrainK)) {
			init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, theta=FALSE, K=constrainK, jumpsize=jumpsize)
		} else {
			init.rate	<- list(values=rep(fit.continuous(ape.tre,orig.dat),length(ape.tre$edge.length)),delta=rep(0,length(ape.tre$edge.length)))
		}
	} else {
		init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, theta=FALSE, K=constrainK, jumpsize=jumpsize )
	}
		
    cur.rates		<- init.rate$values			
	cur.delta.rates	<- init.rate$delta
	cur.root		<- adjust.value(mean(orig.dat), jumpsize)
	cur.vcv			<- update.vcv(ape.tre, cur.rates)
	
	mod.cur = new.bm.lik.fn(cur.rates, cur.root, orig.dat, cur.vcv)
		
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
	prob.rates=1-sum(c(2*prob.mergesplit, prob.root))
	if(prob.rates<0) stop("proposal frequencies must not exceed 1; adjust 'prob.mergesplit' and (or) 'prob.root'")
	proposal.rates<-orig.proposal.rates<-c(2*prob.mergesplit, prob.root, prob.rates)
	prop.cs<-orig.prop.cs<-cumsum(proposal.rates)
	tickerFreq=ceiling(ngen/30)
	tuneFreq.tmp=unique(ceiling(10^(1:ceiling(log(ngen,10)))))
	tuneFreq=tuneFreq.tmp[tuneFreq.tmp<=(0.1*ngen)]
	cur.acceptancerate=0
	
# file handling
	if(summary) parmBase=paste(model, fileBase, "parameters/",sep=".") else parmBase=paste(".", runif(1), model, fileBase, "parameters/",sep="")
	if(!file.exists(parmBase)) dir.create(parmBase)
	parlogger(model=model, init=TRUE, node.des, parmBase=parmBase)
	errorLog=file(paste(parmBase,paste(cur.model, fileBase, "rjTraitErrors.log",sep="."),sep="/"),open='w+')
	generate.error.message(i=NULL, mod.cur=NULL, mod.new=NULL, lnR=NULL, errorLog=errorLog, init=TRUE)
	runLog=file(paste(parmBase,paste(cur.model, fileBase, "rjTrait.log",sep="."),sep="/"),open='w+')
	generate.log(bundled.parms=NULL, cur.model, file=runLog, init=TRUE)

#	orig.internal.only=internal.only

### Begin rjMCMC
    for (i in 1:ngen) {
		
## added 12.16.2010 ## JME v
#		if(i < round(0.1*ngen)) {
#			prop.cs=cumsum(c(0.5, 0.01, 1-(sum(c(0.5, prob.root)))))
#			internal.only=TRUE
#		} else {
#			prop.cs=orig.prop.cs
#			internal.only=orig.internal.only
#		}
## added 12.16.2010 ## JME ^

        lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
		startparms = c(nRateProp, nRateSwapProp, nRcatProp, nRootProp)
		
		if(internal.only) {
				tips.tmp=which(as.numeric(names(cur.delta.rates))<=Ntip(ape.tre))
				if(any(cur.delta.rates[tips.tmp]==1)) stop("Broken internal only sampling")
		}

		
## BM IMPLEMENTATION ##
		while(1) {
			cur.proposal=min(which(runif(1)<prop.cs))
			if (cur.proposal==1 & !constrainK) {									# adjust rate categories
				nr=split.or.merge(cur.delta.rates, cur.rates, ape.tre, node.des, lambdaK, theta=FALSE, internal.only)
				new.rates=nr$new.values
				new.delta.rates=nr$new.delta
				new.root=cur.root
				nRcatProp=nRcatProp+1
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
				break()
			} else if(cur.proposal==2) {											# adjust root
				new.root=adjust.value(cur.root, jumpsize)
				new.rates=cur.rates
				new.delta.rates=cur.delta.rates
				nRootProp=nRootProp+1
				break()
			} else if(cur.proposal==3){													
				if(runif(1)>0.05 & sum(cur.delta.rates)>1) {						# tune single rate
					new.rates=tune.rate(cur.rates, jumpsize) ##
					new.delta.rates=cur.delta.rates 
					new.root=cur.root
					nRateSwapProp=nRateSwapProp+1
					break()
				} else {															# adjust rates
					new.rates=adjust.rate(cur.rates, jumpsize)
					new.delta.rates=cur.delta.rates
					new.root=cur.root
					nRateProp=nRateProp+1
					break()
				} 
			}
		}
		if(any(new.rates!=cur.rates)) new.vcv=update.vcv(ape.tre, new.rates) else new.vcv=cur.vcv
		mod.new=try(new.bm.lik.fn(new.rates, new.root, orig.dat, new.vcv), silent=TRUE)
		
		if(inherits(mod.new, "try-error")) {mod.new=as.list(mod.new); mod.new$lnL=Inf}
		if(!is.infinite(mod.new$lnL)) {
			lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
		} else {
			lnLikelihoodRatio = -Inf
			failures=paste("./failed", model, "parameters/",sep=".")
			if(!file.exists(failures)) dir.create(failures)
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
			mod.cur <- mod.new
			cur.vcv <- new.vcv
			curr.lnL <- mod.new$lnL
			cOK <- determine.accepted.proposal(startparms, endparms, cOK)
		} else {						# deny proposal	
			decision="reject"
			curr.lnL <- mod.cur$lnL 
		}
		
	# potential errors
		if(is.infinite(curr.lnL)) stop("starting point has exceptionally poor likelihood")
		if(r$error) generate.error.message(i, mod.cur, mod.new, lnR, errorLog)
		
	# iteration-specific functions
		if(i%%tickerFreq==0 & summary) {
			if(i==tickerFreq) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
			cat(". ")
		}
		
		if(i%%sample.freq==0) {
			if(colorplots) branchcol.plot(phy,cur.rates,show.tip.label=FALSE)
			bundled.parms=list(gen=i, mrate=exp(mean(log(cur.rates))), cats=sum(cur.delta.rates)+1, root=cur.root, lnL=curr.lnL)
			generate.log(bundled.parms, cur.model, file=runLog)			
			parlogger(model=model, init=FALSE, node.des=node.des, i=i, curr.lnL=curr.lnL, cur.root=cur.root, cur.rates=cur.rates, cur.delta.rates=cur.delta.rates, parmBase=parmBase)
		}
									   
		if(any(tuneFreq%in%i) & adjustable.jumpsize) {
			new.acceptancerate=sum(cOK)/sum(endparms)
			if((new.acceptancerate<cur.acceptancerate) & adjustable.jumpsize) {
				jumpsize=calibrate.jumpsize(phy, data, jumpsizes=c(exp(log(jumpsize)-log(2)), exp(log(jumpsize)+log(2))))
			}
			cur.acceptancerate=new.acceptancerate
		}
	}
	
# End rjMCMC
	close(errorLog)
	close(runLog)
	if(summary) summarize.run(cOK, endparms, cur.model)	else unlink(parmBase, recursive=TRUE)
	return(list(acceptance.rate=sum(cOK)/sum(endparms), jumpsize=jumpsize))
}

rjMCMC.ou<-function (phy, data, ngen=1000, sampleFreq=100, probMergeSplit=0.05, probRoot=0.01, probRate=0.1, probAlpha=0.01, upperAlpha=100, lambdaK=log(2), constrainK=FALSE, jumpsize=NULL, fileBase="result", simplestart=FALSE, switcheroo=TRUE) 
{
#	load("rjMCMC/versions/fakedata.rda"); phy=sims$phy; data=sims[[2]]; ngen=1000; sampleFreq=100; probMergeSplit=0.05; probRoot=0.01; probRate=0.1; probAlpha=0.01; lambdaK=log(2); constrainK=FALSE; jumpsize=2; upperAlpha=100; fileBase="result"; simplestart=FALSE 
	
	model="OU"
	heat=1
	if(is.null(jumpsize)) jumpsize=log(Ntip(phy)/3)
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
	
	cur.rate		<- fit.continuous(ape.tre,orig.dat)
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
				nr=split.or.merge(cur.delta.theta, cur.theta, ape.tre, node.des, lambdaK, theta=TRUE)
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
					} else if(runif(1)<0.5*(probSwap/minorprobs) & length(unique(cur.theta))>1 & switcheroo) {	# swap thetas
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
#	summarize.run(cOK, endparms, cur.model)	
	return(sum(cOK)/sum(endparms))
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

new.bm.lik.fn <- function(rates, root, orig.dat, vcv) { # from LJ HARMON	
# mod 12.02.2010 JM Eastman: using determinant()$modulus rather than det() to stay in log-space
	y=orig.dat
	
	b <- vcv
	w <- rep(root, nrow(b))
	num <- -t(y-w)%*%solve(b)%*%(y-w)/2
	den <- 0.5*(length(y)*log(2*pi) + as.numeric(determinant(b)$modulus))
	list(
		 root=root,
		 lnL = (num-den)[1,1]
		 )	
}

split.or.merge <- function(cur.delta, cur.values, phy, node.des, lambda=lambda, theta=FALSE, internal.only=FALSE) { 
# revert to 11.17.2010 split JM Eastman
#	cur.delta=cur.delta.rates; cur.values=cur.rates; theta=FALSE; lambda=log(2); jumpsize=2
	bb=cur.delta
	vv=cur.values
	names(vv)<-names(bb)<-phy$edge[,2]
	new.bb=choose.one(bb, phy=phy, internal.only)
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
		if(theta) u=split.value(cur.vv=cur.vv, n.desc=n.desc, n.split=n.split, factor=lambda) else u=split.rate(cur.vv=cur.vv, n.desc=n.desc, n.split=n.split, factor=lambda)
		nr.split=u$nr.split
		nr.desc=u$nr.desc
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
# excluded 12.11.2010 JM Eastman
#	cur.delta=cur.delta.rates; cur.values=cur.rates; theta=FALSE; lambda=log(2); jumpsize=2

	vv=cur.values
	rates.sample=unique(vv)[sample(1:length(unique(vv)), 2)]
	new.vv=vv
	new.vv[which(vv==rates.sample[1])]=rates.sample[2]
	new.vv[which(vv==rates.sample[2])]=rates.sample[1]
	
	return(list(new.values=new.vv))
}

calibrate.jumpsize=function(phy, dat, nsteps=100, jumpsizes=NULL) {
#	print("CALIBRATING jumpsize...")

	if(!within.range(nsteps, 100, 1000)) {
		if(nsteps>1000) nsteps=1000 else nsteps=100
	}
	if(is.null(jumpsizes)) jumpsizes=2^(-2:3)
	acceptance.rates=sapply(jumpsizes, function(x) rjMCMC.bm(phy, dat, nsteps, jumpsize=x, summary=FALSE, fileBase="jumpsize", internal.only=TRUE)$acceptance.rate)
	return(sum(jumpsizes*(acceptance.rates/sum(acceptance.rates))))
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
		init.rate=fit.continuous(phy,data)
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

design.vcv.phylo=function(phy, array=TRUE) {
	n=Ntip(phy)
	tips=lapply(phy$edge[,2], function(x) get.descendants.of.node(x,phy,tips=TRUE))
	if(!array) {
		dd=matrix(0,n,n)
		design=lapply(as.numeric(phy$edge[,2]), function(x) {
				  new=dd
				  node=x
				  marker=as.numeric(which(phy$edge[,2]==node))
				  if(node<=n) tt=node else tt=tips[[marker]]
				  for(i in length(tt):1) {
					j=tt[1:i]
					new[tt[i],j]<-new[j,tt[i]]<-1
				  }	  
				  return(new)
				  
		}
		)
	} else {
		design=array(0,dim=c(n,n,length(phy$edge[,2])))
		for(k in 1:length(phy$edge[,2])){
			node=as.numeric(phy$edge[k,2])
			marker=k
			if(node<=n) tt=node else tt=tips[[marker]]
			for(i in length(tt):1) {
				j=tt[1:i]
				design[tt[i],j,k]<-design[j,tt[i],k]<-1
			}	  
		}		 
	} 
	return(design) 
}

update.vcv=function(ape.tre, new.rates) {
	n=ape.tre
	n$edge.length=n$edge.length*new.rates
	vv=vcv.phylo(n)
	return(vv)
}

update.vcv.older=function(ape.tre, design.vcv, new.rates) {
	n=Ntip(ape.tre)
	new.vcv=matrix(0,n,n)
	new.edges=ape.tre$edge.length*new.rates
	for(i in 1:length(new.rates)) {
		new.vcv=new.vcv+new.edges[i]*design.vcv[[i]]
	}
	return(new.vcv)
}

update.vcv.new=function(ape.tre, new.rates, scalar.list) {
#	phy=birthdeath.tree(1,0,taxa=200); d=design.vcv.phylo(phy); s=build.scalar.sets(d); Rprof(); system.time(update.vcv.new(phy, rep(1,length(phy$edge.length)), s)); print(summaryRprof()); system.time(vcv.phylo(phy))
	n=Ntip(ape.tre)
	e=new.rates*ape.tre$edge.length
	ss=scalar.list*e
	entries=apply(ss,3,sum)
	return(matrix(entries,n,n))
}

update.vcv.newest=function(ape.tre, new.rates, scalar.list) {
#	phy=birthdeath.tree(1,0,taxa=200); d=design.vcv.phylo(phy); s=build.scalar.sets(d); Rprof(); system.time(update.vcv.new(phy, rep(1,length(phy$edge.length)), s)); print(summaryRprof()); system.time(vcv.phylo(phy))
	n=Ntip(ape.tre)
	e=new.rates*ape.tre$edge.length
	ss=scalar.list$scalars*e
	new.vcv=matrix(0,n,n)
	dg=which(scalar.list$diags)
	cv=which(scalar.list$attr)
	new.vcv[dg]=apply(scalar.list$scalars[,,dg],2,sum)
	new.vcv[cv]=apply(scalar.list$scalars[,,cv],2,sum)
	return(new.vcv)
}


build.scalar.sets=function(design.vcv) {
	if(!is.array(design.vcv)) stop("Cannot process design matrix: must be array")
	d=design.vcv
	t=dim(design.vcv)[1]
	ss=array(dim=c(1,(2*t-2),(t^2)))
	index=0
	scalar.attributes=c()
	diags=c()
	for(i in 1:t) {
		for(j in 1:t) {
			index=index+1
			ss[,,index]=(as.logical(d[i,j,]))
			if(any(ss[,,index])) scalar.attributes[index]=TRUE else scalar.attributes[index]=FALSE
			if(i==j) diags[index]=TRUE else diags[index]=FALSE
		}
	}
	for(i in 1:index) {
		
	}
	return(list(scalars=ss, attr=scalar.attributes, diags=diags))
}

get.descendants.of.node <- function(node, phy, tips=FALSE) {
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
		storage=as.numeric(sort(unique(unlist(union(desc,storage)))))
	}
	if(tips) return(storage[storage<=Ntip(phy)]) else return(storage)
}

get.desc.of.node <- function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}

get.ancestor.of.node<-function(node, phy) {
	anc=phy$edge[which(phy$edge[,2]==node),1]
	anc
}

get.ancestors.of.node<-function(node, phy) {
	ancestors=c()
	a=0
	while(1) {
		a=a+1
		node<-as.numeric(get.ancestor.of.node(node, phy))
		if(!length(node)) break() else ancestors[a]=node
	}
	return(ancestors)
}

choose.one <- function(cur.delta, phy=NULL, internal.only=FALSE)
# updated 12.16.2010
{
	bb=cur.delta
	if(internal.only) {
		while(1) {
			s=sample(1:length(bb),1)
#			print(bb[s])
			names(bb)=phy$edge[,2]
			if(as.numeric(names(bb[s]))>Ntip(phy)) break()
		}
	} else {
		s=sample(1:length(bb),1)
	}
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

tune.value <- function(values, jumpsize) {
	ss=sample(values, 1)
	ww=which(values==ss)
	nn=adjust.value(ss,jumpsize)
	values[ww]=nn
	return(values)
}

adjust.rate <- function(rate, jumpsize) {
# mod 10.20.2010 JM Eastman
	vv=rate
	v=log(vv)
	rch <- runif(1, min = -jumpsize, max = jumpsize)
	v=exp(v[]+rch)
	v
}

tune.rate <- function(rates, jumpsize) {
	ss=sample(rates, 1)
	ww=which(rates==ss)
	nn=adjust.rate(ss,jumpsize)
	rates[ww]=nn
	return(rates)
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
				 

split.rate=function(cur.vv, n.desc, n.split, factor=log(2)){
	dev=cur.vv-exp(log(cur.vv)+runif(1, -factor, factor))
	nr.desc=cur.vv+dev/n.desc
	nr.split=cur.vv-dev/n.split
#	print((nr.split*n.split+nr.desc*n.desc)/(n.desc+n.split)-cur.vv)
	return(list(nr.desc=nr.desc, nr.split=nr.split))
}

split.value=function(cur.vv, n.desc, n.split, factor=log(2)){
	dev=cur.vv-adjust.value(cur.vv, factor)
	nr.desc=cur.vv + dev/n.desc
	nr.split=cur.vv - dev/n.split
#	print((nr.split*n.split+nr.desc*n.desc)/(n.desc+n.split)-cur.vv)
	return(list(nr.desc=nr.desc, nr.split=nr.split))
}

collect.nodes=function(tips,phy,strict=FALSE){
	if(length(tips)==1) stop("Must supply at least two tips.")
	tt=phy$tip.label
	tt=lapply(tips, function(x) which(phy$tip.label==x))
	nodes=lapply(tt, function(x) get.ancestors.of.node(x, phy))
	all.nodes=nodes[[1]]
	for(i in 2:length(nodes)) {
		all.nodes=intersect(all.nodes,nodes[[i]])
	}
	base.node=max(all.nodes)
	if(strict) {
		u=unlist(nodes)
		return(c(sort(unique(u[u>=base.node])), unlist(tt)))
	} else {
		return(c(base.node,get.descendants.of.node(base.node,phy)))
	}
}

find.relatives=function(tips,phy){
	nn=collect.nodes(tips, phy, strict=FALSE)
	return(phy$tip.label[nn[nn<=Ntip(phy)]])
}

fit.continuous=function(phy, dat){
		n=Ntip(phy)
		ic=try(pic(dat,phy),silent=TRUE)
		if(!inherits(ic, "try-error")) {
			r=mean(ic^2)*((n-1)/n)
			return(r)
		} else {
			return(fitContinuous(phy,dat)$Trait1$beta)
		}
}

## SUMMARIZING MCMC -- PLOTTING FUNCTIONS ##
shifts.plot=function(phy, base.dir, burnin=0, level=0.03, internal.only=FALSE, paint.branches=TRUE, pdf=TRUE, verbose=TRUE, lab="OUT", ...) {
#	phy; base.dir; burnin=0; level=0.03; internal.only=FALSE; paint.branches=TRUE; pdf=TRUE; verbose=TRUE; lab="OUT"
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
		colors.branches.tmp=branchcol.plot(phy, ests, plot=FALSE)
		colors.branches=colors.branches.tmp$col
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
		names(colors.tips)=1:Ntip(phy)
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
#	nn=(Ntip(phy)+1):max(phy$edge[,2])
	NN=phy$edge[,2]
#	ll<-cc<-rr<-rep(0,length(nn))
	ll<-cc<-rr<-rep(0,length(NN))
	if(length(nodes)) {
		ll[NN%in%nodes]=1
		cc[match(names(aa), NN)]=aa
		rr[match(names(colors.nodes), NN)]=colors.nodes
		names(cols)=names(aa)
		names(desc)=paste("node",nodes,sep=".")
	}
	if(length(aa.tips)) {
		tips=as.numeric(names(aa.tips))
		ll[NN%in%tips]=1
		cc[match(names(aa.tips), NN)]=aa.tips
		rr[match(names(colors.tips.tmp), NN)]=colors.tips.tmp
	}
	edgelabels.rj(pch=ifelse(ll==1, 21, NA), cex=2*cc, col=rr, bg=rr, lwd=0.5)
#	ll[nn%in%nodes]=1
#	if(length(nodes)) cc[nn%in%nodes]=aa[match(nn[nn%in%nodes],nodes)]/max(c(aa,tt.cex)) 
#	if(length(nodes)) cc[nn%in%nodes]=aa[match(nn[nn%in%nodes],nodes)]
#	if(!is.null(colors.nodes)) rr[nn%in%nodes]=colors.nodes[order(names(colors.nodes))]
#	nodelabels(pch=ifelse(ll==1, 21, NA), cex=2*cc, col=rr, bg=rr, lwd=0.5)
	
#	tiplabels(pch=ifelse(tt==1, 21, NA), cex=2*tt.cex, col=colors.tips, bg=colors.tips, lwd=0.5)
	legend("topright", title=toupper("direction"), cex=0.5, pt.cex=1, text.col="darkgray", legend = sprintf("%6.1f", rev(c.seq[seq(1,color.length,by=2)])), pch=21, ncol=1, col = "darkgray", pt.bg = rev(ccx[seq(1,color.length,by=2)]), box.lty="blank", border="white")
#	if(length(nodes) | (length(aa.tips) & !internal.only)) {
	if(any(colors.branches!=1)) legend("right", title=toupper("rates"), cex=0.5, pt.cex=1, text.col="darkgray", legend = sprintf("%.6f", colors.branches.tmp$legend.seq), pch=21, ncol=1, col = "darkgray", pt.bg = rev(colors.branches.tmp$legend.col), box.lty="blank", border="white")
	legend("bottomright", title=toupper("frequency"), cex=0.5, pt.cex=2*(seq(0, 1, length=9)), text.col="darkgray", legend = sprintf("%10.3f", seq(0, 1, length=9)), pch=21, ncol=1, col = "darkgray", pt.bg = "white", box.lty="blank", border="white")
#	}
	if(pdf) dev.off()
	
#	if(length(nodes)) {
#		names(cols)=names(aa)
#		names(desc)=paste("node",nodes,sep=".")
#	}
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

edgelabels.rj =function (text, edge, adj = c(0.5, 0.5), frame = "rect", pch = NULL, 
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
    ape:::BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo, pie, 
					 piecol, col, bg, horiz, width, height, ...)
}

shifts.plot.sim=function(phy, base.dir, burnin=0, level=0.03, internal.only=FALSE, paint.branches=TRUE, pdf=TRUE, verbose=TRUE, lab="OUT", curnodes=NA, ...) {
#	phy; base.dir; burnin=0; level=0.03; internal.only=FALSE; paint.branches=TRUE; pdf=TRUE; verbose=TRUE; lab="OUT"
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
		colors.branches.tmp=branchcol.plot(phy, ests, plot=FALSE)
		colors.branches=colors.branches$col
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
		names(colors.tips)=1:Ntip(phy)
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
#	nn=(Ntip(phy)+1):max(phy$edge[,2])
	NN=phy$edge[,2]
#	ll<-cc<-rr<-rep(0,length(nn))
	ll<-cc<-rr<-nn.ref<-rep(0,length(NN))
	if(length(nodes)) {
		ll[NN%in%nodes]=1
		cc[match(names(aa), NN)]=aa
		rr[match(names(colors.nodes), NN)]=colors.nodes
		names(cols)=names(aa)
		names(desc)=paste("node",nodes,sep=".")
	}
	if(length(aa.tips)) {
		tips=as.numeric(names(aa.tips))
		ll[NN%in%tips]=1
		cc[match(names(aa.tips), NN)]=aa.tips
		rr[match(names(colors.tips.tmp), NN)]=colors.tips.tmp
	}
#	nn=(Ntip(phy)+1):max(phy$edge[,2])
#	nn.ref=rep(0,Nnode(phy))
	if(!is.na(curnodes)) nn.ref[NN%in%curnodes]=1
	edgelabels.rj(pch=ifelse(nn.ref==1, 21, NA), cex=3, col=1, bg=NA)
	
	edgelabels.rj(pch=ifelse(ll==1, 21, NA), cex=2*cc, col=rr, bg=rr, lwd=0.5)
#	ll[nn%in%nodes]=1
#	if(length(nodes)) cc[nn%in%nodes]=aa[match(nn[nn%in%nodes],nodes)]/max(c(aa,tt.cex)) 
#	if(length(nodes)) cc[nn%in%nodes]=aa[match(nn[nn%in%nodes],nodes)]
#	if(!is.null(colors.nodes)) rr[nn%in%nodes]=colors.nodes[order(names(colors.nodes))]
#	nodelabels(pch=ifelse(ll==1, 21, NA), cex=2*cc, col=rr, bg=rr, lwd=0.5)
	
#	tiplabels(pch=ifelse(tt==1, 21, NA), cex=2*tt.cex, col=colors.tips, bg=colors.tips, lwd=0.5)
	legend("topright", title=toupper("direction"), cex=0.5, pt.cex=1, text.col="darkgray", legend = sprintf("%6.1f", rev(c.seq[seq(1,color.length,by=2)])), pch=21, ncol=1, col = "darkgray", pt.bg = rev(ccx[seq(1,color.length,by=2)]), box.lty="blank", border="white")
#	if(length(nodes) | (length(aa.tips) & !internal.only)) {
	legend("bottomright", title=toupper("frequency"), cex=0.5, pt.cex=2*(seq(0, 1, length=9)), text.col="darkgray", legend = sprintf("%10.3f", seq(0, 1, length=9)), pch=21, ncol=1, col = "darkgray", pt.bg = "white", box.lty="blank", border="white")
#	}
	if(pdf) dev.off()
	
#	if(length(nodes)) {
#		names(cols)=names(aa)
#		names(desc)=paste("node",nodes,sep=".")
#	}
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

branchcol.plot=function(phy, cur.rates, color.length=8, digits=3, plot=TRUE, legend=TRUE, legend.title="", ...) {
	if(any(cur.rates<=0)) log=FALSE else log=TRUE
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


## A POSTERIORI ANALYSES
compare.rates=function(nodes=list(A=c(NULL,NULL),B=c(NULL,NULL)), phy, posterior.rates, ...){
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

randomization.test=function(obs=obs, exp=exp, iter=10000, verbose=F, two.tailed=F){
	if(verbose && length(obs)<iter) warning("iterations exceeds observed data")
	if(verbose) warning(paste("p-value roughly approximated for ",length(exp), " values\n",sep=""))
	
	O=sample(obs, iter, replace=TRUE)
	E=sample(exp, iter, replace=TRUE)
	
	result=list(c(O-E))
	p=array(NA, iter)
	
	for(i in 1:length(result[[1]])){
		if(result[[1]][i]>0) {
			p[i]=1
		} else {
			if(result[[1]][i]==0) {
				p[i]=sample(c(0,1),1)
			} else {
				p[i]=0
			}
		}
	}
	
	p.g=round(1-(sum(p)/iter), digits=7)
	p.l=1-p.g
	
	if(verbose) res=list(c(O-E),greater=p.g,lesser=p.l) else res=list(greater=p.g,lesser=p.l)
	if(two.tailed) res=min(c(p.g,p.l))*2 else res=res
	return(res)
}

trace.plot=function(obj, col, alpha, lwd=1, hpd=0.95, bars=TRUE, legend.control=list(plot=TRUE, pos=NA, cex=1, pt.cex=1, pch=22, title=""), truncate=list(min=NA, max=NA), xlim=list(min=NA, max=NA), ...){
###############
# obj: a vector or dataframe (where each column gets separately plotted; names of columns are retrieved for legend plotting)
# col: either empty or a vector of colors of length ncol(obj) 
# alpha: translucence of colors (1~opaque; 0~clear)
# lwd: line width for barbells and density outline
# hpd: highest density region of data to be shaded
# bars: whether to plot a bar representing the hdr of the data
# legend.control: can specify changes to any (or none of these); see ?legend
# truncate: trim data *before* hdr is calculated (if either is not NA)
# xlim: must be specified as a list, if used; this is solely for plotting purposes, as the hdr (if still visible) is valid for the untrimmed data
# ...: additional options associated with plot(); see ?plot and ?par
###############
	
	if(!is.data.frame(obj)) {
		if(is.vector(obj)) obj=as.data.frame(obj) else stop("Object must be supplied as a vector or dataframe.")
	}
	control=list(plot=TRUE,pos=NA,cex=1,pt.cex=1,pch=22,title="")
	control[names(legend.control)]<-legend.control
	legend.control=control
	if(missing(col)) col=gray.colors(ncol(obj))
	if(length(col)!=ncol(obj)) col=gray.colors(ncol(obj))
	if(missing(alpha)) alpha=0.8
	line.col=add.transparency(col, min(c(0.95, alpha+alpha*0.25)))
	col=add.transparency(col, alpha)
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
	dd=lapply(1:ncol(obj), function(x) {y=range(obj[,x],na.rm=TRUE); density(obj[,x], from=min(y), to=max(y), na.rm=TRUE, n=1024, kernel="cosine")})
	if(!is.null(hpd)) {
		hh=lapply(1:ncol(obj), function(z) {zz=obj[!is.na(obj[,z]),z]; hdr(zz, hpd=hpd)})
	} else {
		hh=lapply(1:ncol(obj), function(z) {zz=obj[!is.na(obj[,z]),z]; c(min(zz), max(zz))})
	}
	xx=lapply(1:ncol(obj),function(z) {zz=obj[!is.na(obj[,z]),z]; return(c(min(zz[which(zz>=hh[[z]][1])]), max(zz[which(zz<=hh[[z]][2])])))})
	density.yy=lapply(1:length(dd), function(z) {sapply(xx[[z]], function(w) {xs=nearest.pair(w, dd[[z]]$x); infer.y(w, dd[[z]]$x[xs], dd[[z]]$y[xs])})})
	
# swap in calculated HDR points (necessary for log-scale plotting with thin densities)
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
	
	ii=lapply(1:length(dd), function(z) sapply(dd[[z]]$x, function(w) within.range(w, min(hh[[z]]), max(hh[[z]]))))
	
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
	
	plot(x=NULL, xlim=lims.x, ylim=lims.y, ylab="density", bty="n", type="n", ...)
	q=seq(0,min(lims.y)+0.4*min(lims.y),length=ncol(obj)+2)
	q=q[-c(1:2)]
	for(i in 1:length(dd)) {
		dat=dd[[i]]
		index=ii[[i]]
		xs=xx[[i]]
		ys=density.yy[[i]]
		polygon(c(dat$x[index], rev(dat$x[index])), c(dat$y[index],rep(0, length(dat$y[index]))), col=col[i], border=col[i])
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

hdr=function(z, hpd, lim=NULL){
	xx <- sort(c(lim, seq(min(z), max(z), length = 1024)))
	ez <- ecdf(z)
	f <- suppressWarnings(approxfun(ez(xx), xx))
	fit <- suppressWarnings(optimize(function(x) f(x + hpd) - f(x), c(0, 1 - hpd)))
	if (inherits(fit, "try-error") || is.na(fit$objective)) stop("HDR interval failure")
	ci <- fit$min
	f(c(ci, ci + hpd))
}

add.transparency=function (col, alpha) {
    tmp <- col2rgb(col)/255
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
}

within.range=function(x, min, max) {			 
	a=sign(x-min)
	b=sign(x-max)
	if(abs(a+b)==2) return(FALSE) else return(TRUE)
}

nearest.pair=function(val, x) {
	dev=abs(x-val)
	names(dev)=1:length(dev)
	return(sort(as.numeric(names(sort(dev))[1:2])))
}

infer.y=function(xx, x, y) {
	f=lm(y~x)
	p=unname(coef(f))
	yy=xx*p[2]+p[1]
	return(yy)
}
tracer=function(base, rpois.lambda=log(2), Bayes.factor=NULL, burnin=0.25, col.line=1, pdf=TRUE, legend=FALSE,...) {
# base is either data.frame (of rates) or a file name (corresponding to log file)
# Bayes.factor should be a vector of two numbers (corresponding to models differing in the number of rates in the tree)
	d=read.table(base,header=TRUE)
	tag.tmp=strsplit(base,".",fixed=TRUE)[[1]]
	tag=paste(tag.tmp[(tag.tmp != "rjTrait" & tag.tmp != "log")],collapse=".")
	
	# PLOTTING
	if(pdf) pdf(paste(tag,"pdf",sep="."))
		# lnL trace
		if(any(names(d)=="lnL")){
			plot(d$gen, d$lnL, main=tag, ylab="log-likelihood", xlab="generation", bty="n", cex=0.5)
		}
	
		# rates
		if(any(names(d)=="mean.rate")) dmr=as.data.frame(d$mean.rate) 
		if(length(col.line)!=ncol(dmr)) {
			col.line=gray.colors(ncol(dmr)+1)
			col.line=col.line[-length(col.line)]
		}
		trace.plot(dmr, col=col.line, xlab="rate", hpd=NULL, legend=FALSE, main="posterior sample of mean rates", ...)	
		#
	
	if(any(names(d)=="rates")) {
		bb=max(1, round(burnin*nrow(d)))
		d=d[bb:nrow(d),]
		prior=rpois(nrow(d), rpois.lambda)+1
		tt.prior=table(prior)
		tt.post=table(d$rates)
		df=array(dim=c(2,max(c(prior,d$rates))))
		for(pr in 1:length(tt.prior)) {
			df[1,(as.numeric(names(tt.prior)[pr]))]=tt.prior[pr]
		}
		for(pp in 1:length(tt.post)) {
			df[2,(as.numeric(names(tt.post)[pp]))]=tt.post[pp]
		}
		df[is.na(df)]=0
		dimnames(df)[[2]]=seq(1,ncol(df))
		barplot(df,col=c("whitesmoke","black"),xlab="rates",ylab="frequency",beside=TRUE)
		legend(x="topright", c("prior","posterior"), fill=c("whitesmoke","black"), bty="n", cex=0.75)
		if(!is.null(Bayes.factor)) {
			is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
			resample <- function(a=c(0,0), b=vector(), resamples=1000, frequency=0.5, multi=FALSE) {
				if(!multi) {
					xx=lapply(1:resamples, function(x) {
							  rr=sample(b, size=round(frequency*length(b)), replace=TRUE)
							  return(sum(rr==a[1])/sum(rr==a[2]))
							  }
							  )
				} else {
					xx=lapply(1:resamples, function(x) {
							  rr=sample(b, size=round(frequency*length(b)), replace=TRUE)
							  return(sum(rr==1)/sum(rr>1))
							  }
							  )
				}
				return(unlist(xx))
			}
			if(length(Bayes.factor)==2 ) {
				if(Bayes.factor[2]!="multi") {
					BF.tmp=sapply(Bayes.factor, function(x)if(sum(d$rates==x)>0)return(TRUE)else return(FALSE))
					if(all(BF.tmp)) {
						BF=sum(d$rates==Bayes.factor[1])/sum(d$rates==Bayes.factor[2])
						resampled.BF=resample(Bayes.factor, d$rates)
						trace.plot(resampled.BF, main=paste("Bayes factor: ", paste(paste(Bayes.factor, "-", sep=""), collapse=" versus "), "rate model of evolution", sep=""), xlab="Bayes factor", legend=FALSE)
						lines(data.frame(c(BF,BF),c(0,max(density(resampled.BF)$y)), lty=4))
					} else {
						BF=paste("Model(s)", paste(Bayes.factor[which(!BF.tmp)], collapse=" and "), "appear(s) to be missing from the posterior sample.", sep=" ")
						resampled.BF=NA
					}
				} else {
					BF.tmp=c((sum(d$rates==1)>0), (sum(d$rates>1)>0))
					if(all(BF.tmp)) {
						BF=sum(d$rates==1)/sum(d$rates>1)
						resampled.BF=resample(Bayes.factor, d$rates, multi=TRUE)
						trace.plot(resampled.BF, main=paste("Bayes factor: ", paste(paste(Bayes.factor, "-", sep=""), collapse=" versus "), "rate model of evolution", sep=""), xlab="Bayes factor", legend=FALSE)
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

pool.rjMCMCsamples=function(base.dirs, lab=""){
	outdir=paste(lab,"combined.rjMCMC",sep=".")
	if(!file.exists(outdir)) dir.create(outdir)
	rates=intercalate.samples(lapply(base.dirs, function(x) {y=read.table(paste(x, dir(x, pattern="rates.txt"), sep="/")); names(y)=as.numeric(y[1,]); return(y[-1,])}))
	rates=rbind(names(rates), rates)
	write.table(rates, paste(outdir, "rates.txt", sep="/"), row.names=FALSE, col.names=FALSE, quote=FALSE)
	
	shifts=intercalate.samples(lapply(base.dirs, function(x) {y=read.table(paste(x, dir(x, pattern="rate.shifts.txt"), sep="/")); names(y)=as.numeric(y[1,]); return(y[-1,])}))
	shifts=rbind(names(shifts), shifts)
	write.table(shifts, paste(outdir, "rate.shifts.txt", sep="/"), row.names=FALSE, col.names=FALSE, quote=FALSE)

	write.table(intercalate.samples(lapply(base.dirs, function(x) {y=read.table(paste(x, dir(x, pattern="rjTrait.log"), sep="/"), header=TRUE); return(y)})), paste(outdir, "rjTrait.log", sep="/"), quote=FALSE, row.names=FALSE)
	write.table(intercalate.samples(lapply(base.dirs, function(x) {y=read.table(paste(x, dir(x, pattern="summary.txt"), sep="/"), header=TRUE); return(y)})), paste(outdir, "summary.txt", sep="/"), quote=FALSE, row.names=FALSE)
}

intercalate.samples.old=function(list.obj) {
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
	r=nrow(list.obj[[1]])
	c=ncol(list.obj[[1]])
	col.names=names(list.obj[[1]])
#	if(length(unique(unlist(lapply(list.obj, nrow))))!=1) stop("Non-matching object lengths")
	out.array=array(dim=c(r*n, c))
	for(nn in 1:n){
		index=nn
		for(rr in 1:r){
			cur.samples=list.obj[[nn]][rr,]
			if(c>1){
				mm=match(names(cur.samples), col.names)
				out.array[index,]=unlist(c(cur.samples[,mm]))
			} else {
				out.array[index]=cur.samples
			}
			index=index+n
		}
	}
	
	if(c==1) dimnames(out.array)[[2]]=list(col.names) else dimnames(out.array)[[2]]=col.names
	return(as.data.frame(out.array))
}

intercalate.samples=function(list.obj) {
	if(length(list.obj)>2) stop("Nothing to intercalate... too few objects in list.")
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
	names(out.array)=orig.names
	return(out.array)
}

