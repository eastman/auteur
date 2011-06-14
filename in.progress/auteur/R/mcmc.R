#function for hidden Markov sampling from a jump-Brownian process of evolution for continuous traits
#author: JM EASTMAN 2011

mcmc.levy <- function (	phy, dat, SE=0, ngen=1000, sample.freq=100, 
						prob.lambda=0.75, prob.sigsqB=0.075, prob.sigsqJ=0.15, 
						prop.width=NULL, lambdaN=log(2), tuner=0.05,
						lim=list(min=0, max=1000), summary=TRUE, fileBase="result") 
{ 
	
#	require(auteurExtended); tmp=phenogram(rcoal(200), "levy", plot=FALSE); phy=tmp$phy; dat=tmp$dat; prob.lambda=0.75; prob.sigsqB=0.075; prob.sigsqJ=0.15; lambdaN=log(2); prop.width=0.5; summary=TRUE; fileBase="result"; ngen=1000; sample.freq=100; SE=0
	model="LEVY"
	primary.parameter="jumps"
	
#	lnl.storage<-subprop.storage<-jump.storage<-c()

# calibration of proposal width	
	if(is.null(prop.width)) {
		cat("CALIBRATING proposal width...\n")
		adjustable.prop.width=TRUE
		prop.width=calibrate.proposalwidth(phy, dat, nsteps=ngen/1000, model=model, lim=lim, widths=NULL)
	} else {
		if(prop.width<=0) stop("please supply a 'prop.width' larger than 0")
		adjustable.prop.width=FALSE
	}
	
### prepare data for hmmMCMC
	dataList		<- prepare.data(phy, dat, SE)			
	ape.tre			<- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat
	SE				<- dataList$SE
	node.des		<- sapply(ape.tre$edge[,2], function(x) {y=get.descendants.of.node(x, ape.tre, tips=TRUE); if(is.null(y)) return(x) else return(y)})
	names(node.des) <- ape.tre$edge[,2]
	
# parameters
	cur.lambda		<- log(2)
	cur.sigsqB		<- fit.continuous(ape.tre, orig.dat)
	cur.sigsqJ		<- cur.sigsqB^(1/4)
	cur.root		<- adjustvalue(mean(orig.dat), prop.width)
	
	vcv				<- cur.jVCV <- vmat(ape.tre)
	cur.jVCV[]		<- 0

# jump table
	cur.jump.table=c()
	if((nj<-rpois(1, cur.lambda)>0)) {
		for(i in 1:nj) {
			cur.jump.table[nj]=c(cur.jump.table, to<-tree.slide(ape.tre)$node)
			cur.jVCV=updatevcv.jumps(cur.jVCV,from=NULL,to=to,node.des)
		}
	}
	cur.nJ		<- length(cur.jump.table)
	cur.vcv		<- updatevcv.levy(vcv, cur.jVCV, cur.sigsqB, cur.sigsqJ)

# compute lnL for starting state
	mod.cur = bm.lik.fn(cur.root, orig.dat, cur.vcv, SE)
	cur.lnL=mod.cur$lnL
					
# proposal frequencies
	prob.root=1-sum(c(prob.lambda, prob.sigsqB, prob.sigsqJ))
	if(prob.root<0) stop("proposal frequencies must not exceed 1; adjust proposal densities.")
	proposal.rates<-orig.proposal.rates<-c(prob.lambda, prob.sigsqB, prob.sigsqJ, prob.root)
	prop.cs<-orig.prop.cs<-cumsum(proposal.rates)

# proposal counts
	n.props<-n.accept<-rep(0,length(prop.cs))
	names.props<-c("lambda","sigsqB","sigsqJ","root")
	cur.acceptancerate=0
		
# file handling
	if(summary) {
		parmBase=paste(model, fileBase, "parameters/",sep=".")
		if(!file.exists(parmBase)) dir.create(parmBase)
		errorLog=paste(parmBase,paste(model, fileBase, "mcmc.errors.log",sep="."),sep="/")
		runLog=file(paste(parmBase,paste(model, fileBase, "mcmc.log",sep="."),sep="/"),open='w+')
		parms=list(gen=NULL, lnL=NULL, lnL.p=NULL, lnL.h=NULL, 
			   delta=NULL, jumps=NULL, lambda=NULL, sigsq.J=NULL, sigsq.B=NULL, root=NULL)
		parlogger(ape.tre, init=TRUE, primary.parameter, parameters=parms, parmBase=parmBase, runLog, class="hmm")
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
	
### Begin MCMC
    for (i in 1:(ngen+max(tuneFreq))) {
		
        lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
		
		if(cur.nJ<0) stop("Negative cur.nJ")
		
## jump-diffusion IMPLEMENTATION ##
		while(1) {
			cur.proposal=min(which(runif(1)<prop.cs))
			subprop<-tmp<-NULL
			if (cur.proposal==1) {													## adjust JUMP parameters ##
				if(runif(1)>0.75 ) {													# adjust lambda
					nr=tune.rate(rates=cur.lambda, prop.width=prop.width, tuner=tuner, lim=lim)
					new.lambda=nr$values
					new.sigsqB=cur.sigsqB
					new.sigsqJ=cur.sigsqJ
					new.root=cur.root
					new.nJ=cur.nJ
					new.jVCV=cur.jVCV
					new.jump.table=cur.jump.table

					# NEED TO CONSIDER PRIOR ODDS (EXPONENTIAL on lambda with shape lambdaN??)
					lnPriorRatio=dexp(new.lambda, lambdaN, log=TRUE) - dexp(cur.lambda, lambdaN, log=TRUE)
					lnHastingsRatio=nr$lnHastingsRatio
					subprop="lambda"
					break()
				} else if(runif(1)>0.5 & length(cur.jump.table)>0){						# adjust location of a jump
					new.lambda=cur.lambda													# -- no change if jump stays on branch
					new.sigsqB=cur.sigsqB
					new.sigsqJ=cur.sigsqJ
					new.root=cur.root
					new.nJ=cur.nJ
										
					tmp=adjustjump(cur.jump.table, cur.jVCV, add=FALSE, drop=FALSE, swap=TRUE, ape.tre, node.des)
					new.jump.table=tmp$jump.table
					new.jVCV=tmp$jVCV
					
					subprop="push"
					break()
				} else {																# adjust N jumps
					if(cur.nJ==0) {
						if(runif(1)<0.5) {													# -- 0 -> 1
							new.nJ=1
							tmp=adjustjump(cur.jump.table, cur.jVCV, add=TRUE, drop=FALSE, swap=FALSE, ape.tre, node.des)
							new.jump.table=tmp$jump.table
							new.jVCV=tmp$jVCV	
							subprop="increment"
						} else {															# -- 0 -> 0
							new.nJ=cur.nJ
							new.jump.table=cur.jump.table
							new.jVCV=cur.jVCV	
							subprop="null"
						}  
					} else {
						if(runif(1)<0.5){													# -- j -> j + 1
							new.nJ=cur.nJ+1
							tmp=adjustjump(cur.jump.table, cur.jVCV, add=TRUE, drop=FALSE, swap=FALSE, ape.tre, node.des)
							new.jump.table=tmp$jump.table
							new.jVCV=tmp$jVCV							
							subprop="increment"
						} else {															# -- j -> j - 1
							new.nJ=cur.nJ-1
							tmp=adjustjump(cur.jump.table, cur.jVCV, add=FALSE, drop=TRUE, swap=FALSE, ape.tre, node.des)
							new.jump.table=tmp$jump.table
							new.jVCV=tmp$jVCV
							subprop="decrement"
						}
					}
					
					if(subprop=="increment") {
						lnPriorRatio=cur.lambda/(cur.nJ+1)
					} else if(subprop=="decrement") {
						lnPriorRatio=cur.nJ/cur.lambda
					} 
					
					if(is.na(lnPriorRatio)) {
						lnPriorRatio=0
						print(list(cur.lambda, cur.nJ, subprop))
					} 
					
					new.lambda=cur.lambda
					new.sigsqB=cur.sigsqB
					new.sigsqJ=cur.sigsqJ
					new.root=cur.root
										
					break()
				} 
			} else if(cur.proposal==2) {											## adjust sigma-squared BROWNIAN ##
				new.lambda=cur.lambda
				nr=tune.rate(rates=cur.sigsqB, prop.width=prop.width, tuner=tuner, lim=lim)
				new.sigsqB=nr$values
				new.sigsqJ=cur.sigsqJ
				new.root=cur.root
				new.jVCV=cur.jVCV
				new.nJ=cur.nJ
				new.jump.table=cur.jump.table
				subprop="sigsq.Brownian"
				lnHastingsRatio=nr$lnHastingsRatio
				break()
			} else if(cur.proposal==3){												## adjust sigma-squared JUMP ##
				new.lambda=cur.lambda
				new.sigsqB=cur.sigsqB
				nr=tune.rate(rates=cur.sigsqJ, prop.width=prop.width, tuner=tuner, lim=lim)
				new.sigsqJ=nr$values
				new.root=cur.root
				new.jVCV=cur.jVCV
				new.nJ=cur.nJ
				new.jump.table=cur.jump.table
				subprop="sigsq.Jump"
				lnHastingsRatio=nr$lnHastingsRatio
				break()
			} else if(cur.proposal==4) {											## adjust root state ##					
				new.lambda=cur.lambda
				new.sigsqB=cur.sigsqB
				new.sigsqJ=cur.sigsqJ
				new.root=proposal.slidingwindow(cur.root, 2*prop.width, lim=list(min=-Inf, max=Inf))$v	
				new.jVCV=cur.jVCV
				new.nJ=cur.nJ
				new.jump.table=cur.jump.table
				subprop="root.state"
				break()
			} 
		}
		if(any(new.jVCV<0)) stop("Broken updating of jump table.")
		
		n.props[cur.proposal]=n.props[cur.proposal]+1				

	# compute fit of proposed model
		if(subprop!="root.state") new.vcv=updatevcv.levy(vcv, new.jVCV, new.sigsqB, new.sigsqJ) else new.vcv=cur.vcv
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
		if(r$error & summary) generate.error.message(Ntip(ape.tre), i, mod.cur, mod.new, lnLikelihoodRatio, lnPriorRatio, lnHastingsRatio, , , errorLog)
		
		if (runif(1) <= r$r) {			## adopt proposal ##
			new.lambda	-> cur.lambda
			new.sigsqB	-> cur.sigsqB
			new.sigsqJ	-> cur.sigsqJ
			new.root	-> cur.root
			new.jVCV	-> cur.jVCV
			new.nJ		-> cur.nJ
			new.vcv	-> cur.vcv
			mod.new		-> mod.cur
			mod.new$lnL	-> cur.lnL 
			new.jump.table -> cur.jump.table
			
			n.accept[cur.proposal] = n.accept[cur.proposal]+1
		} 
#		cur.lnL -> lnl.storage[i]
#		subprop -> subprop.storage[i]
#		cur.nJ	-> jump.storage[i]
		
		# iteration-specific functions
		if(i%%tickerFreq==0 & summary) {
			if(i==tickerFreq) cat("|",rep(" ",9),toupper("generations complete"),rep(" ",9),"|","\n")
			cat(". ")
		}
		
		if(i%%sample.freq==0 & i>max(tuneFreq) & summary) {			
			parms=list(gen=i-max(tuneFreq), lnL=cur.lnL, lnL.p=lnPriorRatio, lnL.h=lnHastingsRatio, 
					   delta=jumps.edgewise(ape.tre, cur.jump.table), jumps=cur.nJ, lambda=cur.lambda, sigsq.J=cur.sigsqJ, sigsq.B=cur.sigsqB, root=cur.root)
			parlogger(ape.tre, init=FALSE, primary.parameter, parameters=parms, parmBase=parmBase, runLog=runLog, class="hmm")
		}
									   
		if(any(tuneFreq%in%i) & adjustable.prop.width) {
			new.acceptancerate=sum(n.accept)/sum(n.props)
			if((new.acceptancerate<cur.acceptancerate) & adjustable.prop.width) {
				prop.width=calibrate.proposalwidth(ape.tre, orig.dat, nsteps=1000, model, widths=c(exp(log(prop.width)-log(2)), exp(log(prop.width)+log(2))))
			}
			cur.acceptancerate=new.acceptancerate
		}
	} # END FOR-LOOP
	
# End rjMCMC
	
  # clear out prop.width calibration directories if calibration run
	if(summary) {
		close(runLog)
		parlogger(primary.parameter=primary.parameter, parmBase=parmBase, end=TRUE, class="hmm")
		summarize.run(n.accept, n.props, names.props)	
		cleanup.files(parmBase, model, fileBase, ape.tre)
	} 
	
	return(list(acceptance.rate=sum(n.accept)/sum(n.props), prop.width=prop.width))
}

