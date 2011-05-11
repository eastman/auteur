#function for Markov sampling from a range of model complexities under the general process of Brownian motion evolution of continuous traits
#author: LJ HARMON 2009, A HIPP 2009, and JM EASTMAN 2010

rjmcmc.bm <-
function (phy, dat, SE=0, ngen=1000, sample.freq=100, prob.mergesplit=0.05, prob.root=0.01, lambdaK=log(2), constrainK=FALSE, jumpsize=NULL, simplestart=FALSE, internal.only=FALSE, summary=TRUE, fileBase="result") 
{ 
	model="BM"
	heat=1 ## heating not currently implemented
	require(geiger)
	
	if(is.null(jumpsize)) {
		cat("CALIBRATING jumpsize...\n")
		adjustable.jumpsize=TRUE
		jumpsize=calibrate.jumpsize(phy, dat, nsteps=ngen/1000, model)
	} else {
		adjustable.jumpsize=FALSE
	}
	

### prepare data for rjMCMC
	cur.model		<- model
	dataList		<- prepare.data(phy, dat, SE)			
	ape.tre			<- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat
	SE				<- dataList$SE
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))

# initialize parameters
	if(is.numeric(constrainK) & (constrainK > length(ape.tre$edge) | constrainK < 1)) stop("Constraint on rate shifts is nonsensical. Ensure that constrainK is at least 1 and less than the number of available nodes in the tree.")

	if(simplestart | is.numeric(constrainK) | internal.only) {
		if(is.numeric(constrainK)) {
			init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, model, K=constrainK, jumpsize=jumpsize)
		} else {
			init.rate	<- list(values=rep(fit.continuous(ape.tre,orig.dat),length(ape.tre$edge.length)),delta=rep(0,length(ape.tre$edge.length)))
		}
	} else {
		init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, model, K=constrainK, jumpsize=jumpsize )
	}
		
    cur.rates		<- init.rate$values			
	cur.delta.rates	<- init.rate$delta
	cur.root		<- adjustvalue(mean(orig.dat), jumpsize)
	cur.vcv			<- updatevcv(ape.tre, cur.rates)
	
	mod.cur = bm.lik.fn(cur.rates, cur.root, orig.dat, cur.vcv, SE)
		
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
		
# find jumpsize calibration generations
	if(adjustable.jumpsize) {
		tf=c()
		tenths=function(v)floor(0.1*v)
		nn=ngen
		while(1) {
			vv=tenths(nn)
			nn=vv
			if(nn<=1) break	else tf=c(nn,tf)		
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
				nr=splitormerge(cur.delta.rates, cur.rates, ape.tre, node.des, lambdaK, theta=FALSE, internal.only)
				new.rates=nr$new.values
				new.delta.rates=nr$new.delta
				new.root=cur.root
				nRcatProp=nRcatProp+1
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
				break()
			} else if(cur.proposal==2) {											# adjust root
				new.root=adjustvalue(cur.root, jumpsize)
				new.rates=cur.rates
				new.delta.rates=cur.delta.rates
				nRootProp=nRootProp+1
				break()
			} else if(cur.proposal==3){													
				if(runif(1)>0.05 & sum(cur.delta.rates)>1) {						# tune local rate
					new.rates=tune.rate(cur.rates, jumpsize)
					new.delta.rates=cur.delta.rates 
					new.root=cur.root
					nRateSwapProp=nRateSwapProp+1
					break()
				} else {															# adjust rates
					new.rates=adjustrate(cur.rates, jumpsize)
					new.delta.rates=cur.delta.rates
					new.root=cur.root
					nRateProp=nRateProp+1
					break()
				} 
			}
		}
		if(any(new.rates!=cur.rates)) new.vcv=updatevcv(ape.tre, new.rates) else new.vcv=cur.vcv
		mod.new=NULL
		mod.new=try(bm.lik.fn(new.rates, new.root, orig.dat, new.vcv, SE), silent=TRUE)
		
		if(inherits(mod.new, "try-error")) {mod.new=as.list(mod.new); mod.new$lnL=-Inf}
		if(!is.infinite(mod.new$lnL)) {
			lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
		} else {
			mod.new$lnL=-Inf
			lnLikelihoodRatio = -Inf
		}
		
	# compare likelihoods
		endparms = c(nRateProp=nRateProp, nRateSwapProp=nRateSwapProp, nRcatProp=nRcatProp, nRootProp=nRootProp)

		r=assess.lnR((heat * lnLikelihoodRatio + heat * lnPriorRatio + lnHastingsRatio)->lnR)

		# potential errors
		if(is.infinite(mod.cur$lnL)) stop("starting point has exceptionally poor likelihood")
		if(r$error) generate.error.message(Ntip(phy), i, mod.cur, mod.new, lnLikelihoodRatio, lnPriorRatio, lnHastingsRatio, cur.delta.rates, new.delta.rates, errorLog)
		
		if (runif(1) <= r$r) {			## adopt proposal ##
			decision="adopt"
			cur.root <- new.root
			cur.rates <- new.rates
			cur.delta.rates <- new.delta.rates
			mod.cur <- mod.new
			cur.vcv <- new.vcv
			curr.lnL <- mod.new$lnL
			cOK <- determine.accepted.proposal(startparms, endparms, cOK)
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
			bundled.parms=list(gen=i-max(tuneFreq), mrate=exp(mean(log(cur.rates))), cats=sum(cur.delta.rates)+1, root=cur.root, lnL=curr.lnL, rates=cur.rates)
			generate.log(bundled.parms, cur.model, file=runLog)			
			parlogger(model=model, init=FALSE, node.des=node.des, i=i-max(tuneFreq), curr.lnL=curr.lnL, cur.root=cur.root, cur.rates=cur.rates, cur.delta.rates=cur.delta.rates, parmBase=parmBase)
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

