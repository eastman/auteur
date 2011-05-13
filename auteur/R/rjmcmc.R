

#function for Markov sampling from a range of model complexities under the general process of Brownian motion evolution of continuous traits
#author: LJ HARMON 2009, A HIPP 2009, and JM EASTMAN 2010-2011

rjmcmc.bm <- function (	phy, dat, SE=0, ngen=1000, sample.freq=100, 
						prob.mergesplit=0.20, prob.root=0.05, lambdaK=log(2), tuner=0.05, 
						constrainK=FALSE, prop.width=NULL, simplestart=FALSE, 
						summary=TRUE, fileBase="result") 
{ 
	model="BM"
	primary.parameter="rates"

# calibration of proposal width	
	if(is.null(prop.width)) {
		cat("CALIBRATING proposal width...\n")
		adjustable.prop.width=TRUE
		prop.width=calibrate.proposalwidth(phy, dat, nsteps=ngen/1000, model)
	} else {
		if(prop.width<=0) stop("please supply a 'prop.width' larger than 0")
		adjustable.prop.width=FALSE
	}
		
### prepare data for rjMCMC
	dataList		<- prepare.data(phy, dat, SE)			
	ape.tre			<- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat
	SE				<- dataList$SE
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))

# initialize parameters
	if(is.numeric(constrainK) & (constrainK > (length(ape.tre$edge)-1) | constrainK < 1)) stop(paste("Constraint on ",primary.parameter, " is nonsensical. Ensure that constrainK is at least 1 and less than the number of edges in the tree.",sep=""))

	if(simplestart | is.numeric(constrainK)) {
		if(is.numeric(constrainK)) {
			init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, logspace=TRUE, K=constrainK, prop.width=prop.width)
		} else {
			init.rate	<- list(values=rep(fit.continuous(ape.tre,orig.dat),length(ape.tre$edge.length)),delta=rep(0,length(ape.tre$edge.length)))
		}
	} else {
		init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, logspace=TRUE, K=constrainK, prop.width=prop.width )
	}
	
#	if(is.null(prior.rate)) prior.rate	<- 10
		
    cur.rates		<- init.rate$values			
	cur.delta.rates	<- init.rate$delta
	cur.root		<- adjustvalue(mean(orig.dat), prop.width)
	cur.vcv			<- updatevcv(ape.tre, cur.rates)
	
	mod.cur = bm.lik.fn(cur.root, orig.dat, cur.vcv, SE)
	cur.lnL=mod.cur$lnL

			
# proposal frequencies
	prob.rates=1-sum(c(2*prob.mergesplit, prob.root))
	if(prob.rates<0) stop("proposal frequencies must not exceed 1; adjust 'prob.mergesplit' and (or) 'prob.root'")
	proposal.rates<-orig.proposal.rates<-c(2*prob.mergesplit, prob.root, prob.rates)
	prop.cs<-orig.prop.cs<-cumsum(proposal.rates)

# proposal counts
	n.props<-n.accept<-rep(0,length(prop.cs))
	names.props<-c("mergesplit","rootstate","rates")
	
	names.subprops<-c("mergesplit","rootstate","ratetune","moveshift","ratescale")
	n.subprops<-n.subaccept<-rep(0,length(names.subprops))
	names(n.subprops)<-names(n.subaccept)<-names.subprops

	cur.acceptancerate=0
		
# file handling
	if(summary) {
		parmBase=paste(model, fileBase, "parameters/",sep=".")
		if(!file.exists(parmBase)) dir.create(parmBase)
		errorLog=paste(parmBase,paste(model, fileBase, "rjmcmc.errors.log",sep="."),sep="/")
		runLog=file(paste(parmBase,paste(model, fileBase, "rjmcmc.log",sep="."),sep="/"),open='w+')
		parms=list(gen=NULL, lnL=NULL, lnL.p=NULL, lnL.h=NULL, 
			   rates=NULL, delta=NULL, root=NULL)
		parlogger(ape.tre, init=TRUE, primary.parameter, parameters=parms, parmBase=parmBase, runLog)
	}
	
# prop.width calibration points
	if(adjustable.prop.width){
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
	
	tickerFreq=ceiling((ngen+max(tuneFreq))/30)
	
### Begin rjMCMC
    for (i in 1:(ngen+max(tuneFreq))) {
		
        lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
		
## BM IMPLEMENTATION ##
		while(1) {
			cur.proposal=min(which(runif(1)<prop.cs))
			if (cur.proposal==1) {													# adjust rate categories
				nr=splitormerge(cur.delta.rates, cur.rates, ape.tre, node.des, lambdaK, logspace=TRUE, internal.only=FALSE)
				new.rates=nr$new.values
				new.delta.rates=nr$new.delta
				new.root=cur.root
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
				subprop="mergesplit"
				break()
			} else if(cur.proposal==2) {											# adjust root
				new.root=proposal.slidingwindow(cur.root, 2*prop.width, min=NULL, max=NULL)$v								
				new.rates=cur.rates
				new.delta.rates=cur.delta.rates
				subprop="rootstate"
				break()
			} else if(cur.proposal==3){													
				if(runif(1)>0.05 & sum(cur.delta.rates)>1) {						# tune local rate
					if(runif(1)>0.5) {
						nr=tune.rate(rates=cur.rates, prop.width=prop.width, min=0, max=Inf, tuner=tuner)
						new.rates=nr$values
						new.delta.rates=cur.delta.rates 
						new.root=cur.root
						lnHastingsRatio=nr$lnHastingsRatio
#						lnPriorRatio=lnprior.exp(cur.rates, new.rates, prior.rate)
						subprop="ratetune"
						break()						
					} else {
						nr=adjustshift(phy=ape.tre, cur.delta=cur.delta.rates, cur.values=cur.rates)
						new.rates=nr$new.values
						new.delta.rates=nr$new.delta 
						new.root=cur.root
						lnHastingsRatio=nr$lnHastingsRatio
						subprop="moveshift"
						break()						
					}
				} else {															# scale rates
					nr=proposal.multiplier(cur.rates, prop.width)
					new.rates=nr$v
					new.delta.rates=cur.delta.rates
					new.root=cur.root
					lnHastingsRatio=nr$lnHastingsRatio
#					lnPriorRatio=lnprior.exp(cur.rates, new.rates, prior.rate)
					subprop="ratescale"
					break()
				} 
			}
		}
		
		n.props[cur.proposal]=n.props[cur.proposal]+1
		n.subprops[subprop]=n.subprops[subprop]+1				


	# compute fit of proposed model
		if(any(new.rates!=cur.rates)) new.vcv=updatevcv(ape.tre, new.rates) else new.vcv=cur.vcv
		mod.new=NULL
		mod.new=try(bm.lik.fn(new.root, orig.dat, new.vcv, SE), silent=TRUE)
		
		if(inherits(mod.new, "try-error")) {mod.new=as.list(mod.new); mod.new$lnL=-Inf}
		if(!is.infinite(mod.new$lnL)) {
			lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
		} else {
			mod.new$lnL=-Inf
			lnLikelihoodRatio = -Inf
		}
		
	# compare likelihoods
		heat=1
		r=assess.lnR((heat * lnLikelihoodRatio + heat * lnPriorRatio + lnHastingsRatio)->lnR)

		# potential errors
		if(is.infinite(mod.cur$lnL)) stop("starting point has exceptionally poor likelihood")
		if(r$error & summary) generate.error.message(Ntip(ape.tre), i, mod.cur, mod.new, lnLikelihoodRatio, lnPriorRatio, lnHastingsRatio, cur.delta.rates, new.delta.rates, errorLog)
		
		if (runif(1) <= r$r) {			## adopt proposal ##
			cur.root <- new.root
			cur.rates <- new.rates
			cur.delta.rates <- new.delta.rates
			cur.vcv <- new.vcv

			n.accept[cur.proposal] = n.accept[cur.proposal]+1
			n.subaccept[subprop] = n.subaccept[subprop]+1
			mod.cur <- mod.new
			cur.lnL <- mod.new$lnL
		} 
		
		# iteration-specific functions
		if(i%%tickerFreq==0 & summary) {
			if(i==tickerFreq) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
			cat(". ")
		}
		
		if(i%%sample.freq==0 & i>max(tuneFreq) & summary) {			
			parms=list(gen=i-max(tuneFreq), lnL=cur.lnL, lnL.p=lnPriorRatio, lnL.h=lnHastingsRatio, 
					   rates=cur.rates, delta=cur.delta.rates, root=cur.root)
			parlogger(ape.tre, init=FALSE, primary.parameter, parameters=parms, parmBase=parmBase, runLog=runLog)
		}
									   
		if(any(tuneFreq%in%i) & adjustable.prop.width) {
			new.acceptancerate=sum(n.accept)/sum(n.props)
			if((new.acceptancerate<cur.acceptancerate) & adjustable.prop.width) {
				prop.width=calibrate.proposalwidth(ape.tre, orig.dat, nsteps=1000, model, widths=c(exp(log(prop.width)-log(2)), exp(log(prop.width)+log(2))))
			}
			cur.acceptancerate=new.acceptancerate
		}
	}
	
# End rjMCMC
	
  # clear out prop.width calibration directories if calibration run
	if(summary) {
		close(runLog)
		parlogger(primary.parameter=primary.parameter, parmBase=parmBase, end=TRUE)
		summarize.run(n.subaccept, n.subprops, names.subprops)	
		cleanup.files(parmBase, model, fileBase, ape.tre)
	} 
	
	return(list(acceptance.rate=sum(n.accept)/sum(n.props), prop.width=prop.width))
}

