#function for Markov sampling from a range of model complexities under the general process of Brownian motion evolution of continuous traits
#author: JM EASTMAN 2011

rjmcmc.ou <- function(	phy, dat, SE=0, ngen=1000, sample.freq=100, 
						prob.mergesplit=0.1, prob.theta=0.4, prob.alpha=0.2, upper.alpha=100, 
						lambdaK=log(2), constrainK=FALSE, prop.width=NULL, simplestart=FALSE, 
						internal.only=FALSE, summary=TRUE, fileBase="result") 
{ 
	model="OU"
	primary.parameter="optima"
	
# calibration of proposal width	
	if(is.null(prop.width)) {
		cat("CALIBRATING proposal width...\n")
		prop.width=calibrate.proposalwidth(phy, dat, nsteps=ngen/1000, model)
		adjustable.prop.width=TRUE		
	} else {
		if(prop.width<=0) stop("please supply a 'prop.width' larger than 0")
		adjustable.prop.width=FALSE
	}
	
# parameter bounds 
	abs.low.alpha=1e-13
	if(length(upper.alpha)==2 & upper.alpha[1]>abs.low.alpha) {
		alpha.bounds=upper.alpha
	} else if(!is.null(upper.alpha) & is.numeric(upper.alpha)){
		alpha.bounds=c(abs.low.alpha, upper.alpha)
	} else {
		alpha.bounds=c(abs.low.alpha, 100)
	}
	
### prepare data for rjMCMC
	dataList		<- prepare.data(phy, dat, SE)			
	ape.tre			<- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat
	SE				<- dataList$SE
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))
	ancestors		<- lapply(1:Ntip(ape.tre), function(x) c(x,rev(sort(get.ancestors.of.node(x,ape.tre)))))
	epochs			<- get.epochs(ape.tre)
	phyobject		<- get.times(ape.tre)
	vcvSE			<- vmat(ape.tre); if(any(SE!=0)) diag(vcvSE)<-diag(vcvSE)+SE^2
	
# initialize parameters
	if(is.numeric(constrainK) & (constrainK > length(ape.tre$edge) | constrainK < 1)) stop(paste("Constraint on ",primary.parameter, " is nonsensical. Ensure that constrainK is at least 1 and less than the number of available nodes in the tree.",sep=""))
	
	if(simplestart | is.numeric(constrainK) | internal.only) {
		if(is.numeric(constrainK) & !internal.only) {
			init.theta	<- generate.starting.point(orig.dat, ape.tre, node.des, logspace=FALSE, K=constrainK, prop.width=prop.width)
		} else {
			init.theta	<- list(values=rep(mean(orig.dat),length(ape.tre$edge.length)),delta=rep(0,length(ape.tre$edge.length)))
		}
	} else {
		init.theta		<- generate.starting.point(orig.dat, ape.tre, node.des, logspace=FALSE, K=constrainK, prop.width=prop.width)
	}
	
	init.hansen		<- fit.hansen(ape.tre,orig.dat,alpha.bounds)

	cur.theta		<- init.theta$values 
	cur.delta.theta	<- init.theta$delta
	cur.alpha		<- init.hansen$alpha
	cur.rate		<- init.hansen$sigmasq
	
# attempt to ensure valid starting point
	while(1) {
		mod.cur = try(ou.lik.fn(orig.dat, ape.tre, epochs, phyobject, ancestors, vcvSE, cur.rate, cur.alpha, c(-999, cur.theta)),silent=TRUE)
		if(inherits(mod.cur, "try-error")){
			if(!summary)stop()
			cur.alpha=adjustrate(cur.alpha, prop.width)
			cur.rate=adjustrate(cur.rate, prop.width)
		} else {
			cur.lnL=mod.cur$lnL
			break()
		}
	}	
	
# proposal frequencies	
	prob.rate=1-sum(c(2*prob.mergesplit, prob.theta, prob.alpha))
	if(prob.rate<0) stop("proposal frequencies must not exceed 1; adjust proposal frequencies")
	proposal.rates<-orig.proposal.rates<-c(2*prob.mergesplit, prob.theta, prob.alpha, prob.rate)
	prop.cs<-orig.prop.cs<-cumsum(proposal.rates)
	
# proposal counts
	n.props<-n.accept<-rep(0,length(prop.cs))
	names.props<-c("mergesplit","optima","alpha","rate")
	cur.acceptancerate=0
	
# file handling
	if(summary) {
		parmBase=paste(model, fileBase, "parameters/",sep=".")
		if(!file.exists(parmBase)) dir.create(parmBase)
		errorLog=paste(parmBase,paste(model, fileBase, "rjmcmc.errors.log",sep="."),sep="/")
		runLog=file(paste(parmBase,paste(model, fileBase, "rjmcmc.log",sep="."),sep="/"),open='w+')
		parms=list(gen=NULL, lnL=NULL, lnL.p=NULL, lnL.h=NULL, 
			   optima=NULL, delta=NULL, alpha=NULL, rate=NULL)
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
		
		if(internal.only) {
			tips.tmp=which(as.numeric(names(cur.delta.theta))<=Ntip(ape.tre))
			if(any(cur.delta.theta[tips.tmp]==1)) stop("Broken internal only sampling")
		}
				
## OU IMPLEMENTATION ##
		while(1) {
			cur.proposal=min(which(runif(1)<prop.cs))
			if (cur.proposal==1 & !constrainK) {									# adjust theta categories
				nr=splitormerge(cur.delta.theta, cur.theta, ape.tre, node.des, lambdaK, logspace=FALSE, internal.only)
				new.theta=nr$new.values ##
				new.delta.theta=nr$new.delta ##
				new.alpha=cur.alpha
				new.rate=cur.rate 
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
				break()
			} else if(cur.proposal==2){													
				if(runif(1)>0.05 & sum(cur.delta.theta)>1) {						# tune local theta
					new.theta=tune.value(cur.theta, prop.width) ##
					new.delta.theta=cur.delta.theta 
					new.alpha=cur.alpha									
					new.rate=cur.rate 					
					break()
				} else {															# adjust theta
					new.theta=adjustvalue(cur.theta, prop.width) ##
					new.delta.theta=cur.delta.theta
					new.alpha=cur.alpha									
					new.rate=cur.rate 
					break()
				} 
			}  else if(cur.proposal==3) {											# adjust alpha
				while(1) {
					new.alpha=adjustrate(cur.alpha, prop.width) ##	
					if(withinrange(new.alpha, min(alpha.bounds), max(alpha.bounds))) break()
				}
				new.theta=cur.theta
				new.delta.theta=cur.delta.theta
				new.rate=cur.rate 
				break()
			} else if(cur.proposal==4) {											# adjust rate
				new.theta=cur.theta
				new.delta.theta=cur.delta.theta
				new.alpha=cur.alpha									
				new.rate=adjustrate(cur.rate,prop.width) ##
				break()
			}
		}
				
		n.props[cur.proposal]=n.props[cur.proposal]+1				
		
# compute fit of proposed model
		mod.new=NULL
		mod.new=try(ou.lik.fn(orig.dat, ape.tre, epochs, phyobject, ancestors, vcvSE, new.rate, new.alpha, c(-999, new.theta)),silent=TRUE)
		
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
		if(r$error & summary) generate.error.message(Ntip(ape.tre), i, mod.cur, mod.new, lnLikelihoodRatio, lnPriorRatio, lnHastingsRatio, cur.delta.theta, new.delta.theta, errorLog)
		
		if (runif(1) <= r$r) {			## adopt proposal ##
			cur.theta <- new.theta
			cur.delta.theta <- new.delta.theta
			cur.alpha <- new.alpha
			cur.rate <- new.rate
			
			n.accept[cur.proposal] = n.accept[cur.proposal]+1
			
			mod.cur <- mod.new
			cur.lnL <- mod.new$lnL
		} 
		
# iteration-specific functions
		if(i%%tickerFreq==0 & summary) {
			if(i==tickerFreq) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
			cat(". ")
		}
		
		if(i%%sample.freq==0 & i>max(tuneFreq) & summary) {

			## DIAGNOSTICS ##
#			cc=cur.theta
#			names(cc)=phy$edge[,2]
#			if(sum(cc)>=1){branchcol.plot(ape.tre, cc, color.length=sum(cur.delta.theta)+1, log=FALSE)}

			parms=list(gen=i-max(tuneFreq), lnL=cur.lnL, lnL.p=lnPriorRatio, lnL.h=lnHastingsRatio, 
					   optima=cur.theta, delta=cur.delta.theta, alpha=cur.alpha, rate=cur.rate)
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
		summarize.run(n.accept, n.props, names.props)	
		cleanup.files(parmBase, model, fileBase, ape.tre)
	} 
	return(list(acceptance.rate=sum(n.accept)/sum(n.props), prop.width=prop.width))
}



#function for Markov sampling from a range of model complexities under the general process of Brownian motion evolution of continuous traits
#author: LJ HARMON 2009, A HIPP 2009, and JM EASTMAN 2010

rjmcmc.bm <- function (	phy, dat, SE=0, ngen=1000, sample.freq=100, 
						prob.mergesplit=0.05, prob.root=0.05, 
						lambdaK=log(2), constrainK=FALSE, prop.width=NULL, simplestart=FALSE, 
						internal.only=FALSE, summary=TRUE, fileBase="result") 
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
	if(is.numeric(constrainK) & (constrainK > length(ape.tre$edge) | constrainK < 1)) stop(paste("Constraint on ",primary.parameter, " is nonsensical. Ensure that constrainK is at least 1 and less than the number of available nodes in the tree.",sep=""))

	if(simplestart | is.numeric(constrainK) | internal.only) {
		if(is.numeric(constrainK) & !internal.only) {
			init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, logspace=TRUE, K=constrainK, prop.width=prop.width)
		} else {
			init.rate	<- list(values=rep(fit.continuous(ape.tre,orig.dat),length(ape.tre$edge.length)),delta=rep(0,length(ape.tre$edge.length)))
		}
	} else {
		init.rate	<- generate.starting.point(orig.dat, ape.tre, node.des, logspace=TRUE, K=constrainK, prop.width=prop.width )
	}
	
	prior.rate		<- fit.continuous(ape.tre,orig.dat)
		
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
		
		if(internal.only) {
				tips.tmp=which(as.numeric(names(cur.delta.rates))<=Ntip(ape.tre))
				if(any(cur.delta.rates[tips.tmp]==1)) stop("Broken internal only sampling")
		}

		
## BM IMPLEMENTATION ##
		while(1) {
			cur.proposal=min(which(runif(1)<prop.cs))
			if (cur.proposal==1 & !constrainK) {									# adjust rate categories
				nr=splitormerge(cur.delta.rates, cur.rates, ape.tre, node.des, lambdaK, logspace=TRUE, internal.only)
				new.rates=nr$new.values
				new.delta.rates=nr$new.delta
				new.root=cur.root
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio + priorratio.exp(cur.rates, cur.delta, new.rates, new.delta, prior.rate)	## MOD 05.06.2011
				break()
			} else if(cur.proposal==2) {											# adjust root
				new.root=proposal.slidingwindow(cur.root, prop.width, min=NULL, max=NULL)$v								## MOD 05.06.2011
				new.rates=cur.rates
				new.delta.rates=cur.delta.rates
				break()
			} else if(cur.proposal==3){													
				if(runif(1)>0.05 & sum(cur.delta.rates)>1) {						# tune local rate
					new.rates=tune.rate(cur.rates, prop.width)
					new.delta.rates=cur.delta.rates 
					new.root=cur.root
					break()
				} else {															# adjust rates
					new.rates=adjustrate(cur.rates, prop.width)
					new.delta.rates=cur.delta.rates
					new.root=cur.root
					break()
				} 
			}
		}
		
		n.props[cur.proposal]=n.props[cur.proposal]+1				

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
		summarize.run(n.accept, n.props, names.props)	
		cleanup.files(parmBase, model, fileBase, ape.tre)
	} 
	
	return(list(acceptance.rate=sum(n.accept)/sum(n.props), prop.width=prop.width))
}

