#function for Markov sampling from a range of model complexities under the general process of Brownian motion evolution of continuous traits
#author: LJ HARMON 2009, A HIPP 2009, and JM EASTMAN 2010

rjmcmc.ou <-
function (phy, dat, SE=0, ngen=1000, sample.freq=100, prob.mergesplit=0.05, prob.root=0.01, prob.rate=0.1, prob.alpha=0.01, upper.alpha=100, lambdaK=log(2), constrainK=FALSE, jumpsize=NULL, simplestart=FALSE, internal.only=FALSE, summary=TRUE, fileBase="result") 
{ 
	model="OU"
	heat=1 ## heating not currently implemented
	require(geiger)
	
	if(is.null(jumpsize)) {
		cat("CALIBRATING jumpsize...\n")
		adjustable.jumpsize=TRUE
		jumpsize=calibrate.jumpsize(phy, dat, nsteps=ngen/1000, model)
	} else {
		adjustable.jumpsize=FALSE
	}
	
	abs.low.alpha=1e-13
	if(length(upper.alpha)==2 & upper.alpha[1]>abs.low.alpha) {
		alpha.bounds=upper.alpha
	} else if(!is.null(upper.alpha) & is.numeric(upper.alpha)){
		alpha.bounds=c(abs.low.alpha, upper.alpha)
	} else {
		alpha.bounds=c(abs.low.alpha, 100)
	}
	
### prepare data for rjMCMC
	cur.model		<- model
	dataList		<- prepare.ouchdata(phy, dat, SE)			
	ape.tre			<- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat
	SE				<- dataList$SE
	ouch.tre		<- dataList$ouch.tre
	ouch.dat		<- dataList$ouch.dat
	ape.to.ouch.order <- match(as.numeric(names(dataList$ape.ordered[-1])), as.numeric(ape.tre$edge[,2]))
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))
	ancestors		<- sapply(1:Ntip(phy), function(x) get.ancestors.of.node(x,phy))
	
# initialize parameters
	if(is.numeric(constrainK) & (constrainK > length(ape.tre$edge) | constrainK < 1)) stop("Constraint on rate shifts is nonsensical. Ensure that constrainK is at least 1 and less than the number of available nodes in the tree.")

	if(simplestart | is.numeric(constrainK) | internal.only) {
		if(is.numeric(constrainK)) {
			init.theta	<- generate.starting.point(orig.dat, ape.tre, node.des, model, K=constrainK, jumpsize=jumpsize)
		} else {
			init.theta	<- list(values=rep(mean(orig.dat),length(ape.tre$edge.length)),delta=rep(0,length(ape.tre$edge.length)))
		}
	} else {
		init.theta		<- generate.starting.point(orig.dat, ape.tre, node.des, model, K=constrainK, jumpsize=jumpsize )
	}
		
    cur.rate		<- fit.continuous(ape.tre, orig.dat)
	cur.alpha		<- exp(runif(1,log(alpha.bounds)))
	cur.theta		<- init.theta$values			
	cur.delta.theta	<- init.theta$delta
	cur.root.theta	<- assignroot(cur.theta,ape.tre)
	cur.regimes		<- c(cur.root.theta,cur.theta) 
	
	mod.cur = ou.lik.fn(cur.rate, cur.alpha, cur.regimes, cur.theta, ouch.dat, ouch.tre, ape.to.ouch.order, SE)	
		
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
	
	prob.theta=1-sum(c(2*prob.mergesplit, prob.root, prob.rate, prob.alpha))
	if(prob.theta<0) stop("proposal frequencies must not exceed 1; adjust proposal frequencies")
	proposal.rates<-orig.proposal.rates<-c(2*prob.mergesplit, prob.root, prob.rate, prob.alpha, prob.theta)
	prop.cs<-orig.prop.cs<-cumsum(proposal.rates)
		
# find jumpsize calibration generations
	if(adjustable.jumpsize) {
		tf=c()
		tenths=function(v)floor(0.1*v)
		nn=ngen
		while(1) {
			vv=tenths(nn)
			nn=vv
			if(nn<=1) break() else tf=c(nn,tf)		
		}
		tuneFreq=tf
	} else {
		tuneFreq=0
	}
	
	cur.acceptancerate=0
	
	tickerFreq=ceiling((ngen+max(tuneFreq))/30)
	
# file handling
	if(summary) parmBase=paste(model, fileBase, "parameters/",sep=".") else parmBase=paste(".", runif(1), model, fileBase, "parameters/",sep="")
	if(!file.exists(parmBase)) dir.create(parmBase)
	parlogger(model=model, init=TRUE, node.des, parmBase=parmBase)
	errorLog=paste(parmBase,paste(cur.model, fileBase, "rjmcmc.errors.log",sep="."),sep="/")
	runLog=file(paste(parmBase,paste(cur.model, fileBase, "rjmcmc.log",sep="."),sep="/"),open='w+')
	generate.log(bundled.parms=NULL, cur.model, file=runLog, init=TRUE)

### Begin rjMCMC
    for (i in 1:(ngen+max(tuneFreq))) {
		branchcol.plot(ape.tre, cur.theta)
		
        lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
		startparms = c(nRateProp, nRootProp, nAlphaProp, nTcatProp, nThetaProp)
		
		if(internal.only) {
				tips.tmp=which(as.numeric(names(cur.delta.theta))<=Ntip(ape.tre))
				if(any(cur.delta.theta[tips.tmp]==1)) stop("Broken internal only sampling")
		}

		
## OU IMPLEMENTATION ##
		while(1) {
			cur.proposal=min(which(runif(1)<prop.cs))
			if (cur.proposal==1 & !constrainK) {									# adjust theta categories
				nr=splitormerge(cur.delta.theta, cur.theta, ape.tre, node.des, lambdaK, theta=TRUE)
				new.alpha=cur.alpha
				new.theta=nr$new.values ##
				new.delta.theta=nr$new.delta ##
				new.root.theta<-nrt<-assignroot(new.theta, ape.tre) ## 
				new.regimes=c(nrt,new.theta) ##
				new.rate=cur.rate 
				nTcatProp=nTcatProp+1
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
				break()
			} else if(cur.proposal==2) {											# adjust root
				new.alpha=cur.alpha									
				new.root.theta=assignroot(cur.theta, ape.tre) ##
				new.theta=cur.theta
				new.delta.theta=cur.delta.theta
				new.regimes=cur.regimes 
				new.rate=cur.rate 
				nRootProp=nRootProp+1					
				break()
			} else if(cur.proposal==3) {											# adjust rate
				new.theta=cur.theta
				new.delta.theta=cur.delta.theta
				new.root.theta=cur.root.theta
				new.alpha=cur.alpha									
				new.regimes=cur.regimes 
				new.rate=adjustrate(cur.rate,jumpsize) ##
				nRateProp=nRateProp+1
				break()
			} else if(cur.proposal==4) {											# adjust alpha
				new.theta=cur.theta
				new.delta.theta=cur.delta.theta
				new.root.theta=cur.root.theta
				new.alpha=adjustrate(cur.alpha, jumpsize) ##								
				new.regimes=cur.regimes 
				new.rate=cur.rate 
				nAlphaProp=nAlphaProp+1
				break()
			} else if(cur.proposal==5){													
				if(runif(1)>0.05 & sum(cur.delta.theta)>1) {						# tune local theta
					new.theta=tune.value(cur.theta, jumpsize) ##
					new.delta.theta=cur.delta.theta 
					new.root.theta=cur.root.theta
					new.alpha=cur.alpha									
					new.regimes=cur.regimes 
					new.rate=cur.rate 					
					nThetaProp=nThetaProp+1
					break()
				} else {															# adjust theta
					new.theta=adjustvalue(cur.theta, jumpsize) ##
					new.delta.theta=cur.delta.theta
					new.root.theta=cur.root.theta
					new.alpha=cur.alpha									
					new.regimes=cur.regimes 
					new.rate=cur.rate 
					nThetaProp=nThetaProp+1
					break()
				} 
			}
		}
		mod.new=NULL
		mod.new=try(ou.lik.fn(new.rate, new.alpha, new.regimes, new.theta, ouch.dat, ouch.tre, ape.to.ouch.order, SE), silent=TRUE)
		
		if(inherits(mod.new, "try-error")) {mod.new=as.list(mod.new); mod.new$lnL=-Inf}
		if(!is.infinite(mod.new$lnL)) {
			lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
		} else {
			mod.new$lnL=-Inf
			lnLikelihoodRatio = -Inf
		}
		
	# compare likelihoods
		endparms = c(nRateProp=nRateProp, nRootProp=nRootProp, nAlphaProp=nAlphaProp, nTcatProp=nTcatProp, nThetaProp=nThetaProp)

		r=assess.lnR((heat * lnLikelihoodRatio + heat * lnPriorRatio + lnHastingsRatio)->lnR)

		# potential errors
		if(is.infinite(mod.cur$lnL)) stop("starting point has exceptionally poor likelihood")
		if(r$error) generate.error.message(Ntip(phy), i, mod.cur, mod.new, lnLikelihoodRatio, lnPriorRatio, lnHastingsRatio, cur.delta.theta, new.delta.theta, errorLog)
		
		if (runif(1) <= r$r) {			## adopt proposal ##
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
		} else {						## deny proposal ##	
			decision="reject" 
			curr.lnL <- mod.cur$lnL 
		}
		
		# iteration-specific functions
		if(i%%tickerFreq==0 & summary) {
			if(i==tickerFreq) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
			cat(". ")
		}
		
		if(i%%sample.freq==0 & i>max(tuneFreq)) {
			bundled.parms=list(gen=i-max(tuneFreq), rate=cur.rate, regimes=sum(cur.delta.theta)+1, root=cur.root.theta, theta=range(cur.theta), alpha=cur.alpha, lnL=curr.lnL)
			generate.log(bundled.parms, cur.model, file=runLog)			
			parlogger(model=model, init=FALSE, node.des=node.des, i=i-max(tuneFreq), curr.lnL=curr.lnL, cur.root=cur.root.theta, cur.rates=cur.rate, cur.alpha=cur.alpha, cur.theta=cur.theta, cur.delta.theta=cur.delta.theta, parmBase=parmBase)
		}
									   
		if(any(tuneFreq%in%i) & adjustable.jumpsize) {
			new.acceptancerate=sum(cOK)/sum(endparms)
			if((new.acceptancerate<cur.acceptancerate) & adjustable.jumpsize) {
				jumpsize=calibrate.jumpsize(phy, dat, nsteps=1000, model, jumpsizes=c(exp(log(jumpsize)-log(2)), exp(log(jumpsize)+log(2))))
			}
			cur.acceptancerate=new.acceptancerate
		}
	}
	
# End rjMCMC
	close(runLog)
	
  # clear out jumpsize calibration directories if calibration run
	if(summary) {
		summarize.run(cOK, endparms, cur.model)	
		cleanup.files(parmBase, cur.model, fileBase)
	} else {
		unlink(parmBase, recursive=TRUE)
	}
	return(list(acceptance.rate=sum(cOK)/sum(endparms), jumpsize=jumpsize))
}

