## Written by LJ HARMON, AH: 2008
## Updated by JM EASTMAN: 09.2010

rjMCMC.trait<-function (phy, data, ngen = 1000, model=c("OU","BM"), probAlpha=0.05, probRegimes=0.10, probMergeSplit = 0.05, heat=1, fileBase = "result") 
{
	if(length(model)>1)stop(paste("Please specify either ", sQuote("OU"), " or ", sQuote("BM")))
	
	require(ouch)
	require(geiger)
	
### prepare objects for rjMCMC
	dataList <- prepare.ouchdata(phy, data)					# check tree and data; convert into S4
	ape.tre <- dataList$ape.tre								# S3 phylogeny
	orig.dat <- dataList$orig.dat							# S3 data
	ouch.tre <- dataList$ouch.tre							# S4 OUCH phylogeny
	ouch.dat <- dataList$ouch.dat							# S4 OUCH data
	nn <- length(ape.tre$edge.length)						# max number of rates (one per branch)
	currRates.c <- numeric(nn)								# initialize array for recording which branch belongs to which rate category
    currRates.c[] <- 1										# assign one rate to all branches initially
    currRates <- numeric(nn)								# initialize array for rate estimates branchwise
    currRates[] <- 0.001									# initialize BMMrate to be small
	currAlpha <- 0.001										# initialize alpha
	currRegimes <- numeric(nrow(ouch.dat))					# initial set of selective regimes
	currRegimes[] <- 1										# initialize regime to be global
	currModel = model										# initialize family of models explored
    nRateProp = 0											# proposal count: update rate parameter(s)
	nRcatProp = 0											# proposal count: update rate categories
	nAlphaProp = 0											# proposal count: update alpha parameter
	nRegimesProp = 0										# proposal count: update selective regimes
	nRateOK = 0												# proposal count: successful rate changes
	nRcatOK = 0												# proposal count: successful rate categories updates
	nAlphaOK = 0											# proposal count: successful alpha updates
	nRegimesOK = 0											# proposal count: successful regime updates
	cOK=c(nRateOK, nRcatOK, nAlphaOK, nRegimesOK)			# 
	l <- array(dim=ngen)									# number of MCMC generations
    l[] <- NA												# initialize array for storing lnLs
	r.nrate<-r.alpha<-r.nregimes<-l
	
	
### Begin rjMCMC
#	cat("\trates\tregimes\t\talpha\t\tlnL\t\t\tmodel")
    for (i in 1:ngen) {
        lnLikelihoodRatio <- lnProposalRatio <- lnPriorRatio <- 0
		startparms=list(currRates, currRates.c, currAlpha, currRegimes)
		model[i]=currModel

## OU IMPLEMENTATION ## thetas currently optimized by ouch:::glssoln(); selective regimes are choosable (see Beta)
		if(currModel=="OU") {
			if (runif(1) < (2 * probMergeSplit)) {
				adj.regimes <- adj.rates <- FALSE
				if(runif(1) < 0.5) adj.regimes=TRUE else adj.rates=TRUE
				if(adj.regimes) {
					nr=split.or.merge(currRegimes, currRegimes, decide.s.or.m(currRegimes)->s.m)
					newRegimes=nr$new.categories
					newRates=currRates
					newRates.c=currRates.c
					nRegimesProp=nRegimesProp+1
				} else { # adjust rate categories
					nr=split.or.merge(currRates, currRates.c, decide.s.or.m(currRates.c)->s.m)
					newRates=nr$new.values
					newRates.c=nr$new.categories
					newRegimes=currRegimes
					nRcatProp=nRcatProp+1
				}
				newAlpha=currAlpha
				lnProposalRatio=nr$lnl.prop
				lnPriorRatio=nr$lnl.prior
			} else { # neither split nor merge
				if(runif(1)<probAlpha) { # adjust alpha
					newAlpha=adjust.rates(currAlpha)
					newRates=currRates
					newRates.c=currRates.c
					nAlphaProp=nAlphaProp+1
				} else { # adjust rates
					newRates=adjust.rates(currRates)
					newRates.c=fix.categories(newRates)
					newAlpha=currAlpha
					nRateProp=nRateProp+1
				}
				newRegimes=currRegimes
			}			
			modCurr=lnl.OU(currRates, currAlpha, currRegimes, ape.tre, ouch.dat)
			modNew=lnl.OU(newRates, newAlpha, newRegimes, ape.tre, ouch.dat)
			
		} else if(currModel=="BM"){
			if (runif(1) < (2 * probMergeSplit)) {
				nr=split.or.merge(currRates, currRates.c, decide.s.or.m(currRates.c)->s.m)
				newRates=nr$new.values
				newRates.c=nr$new.categories
				nRcatProp=nRcatProp+1
				lnProposalRatio=nr$lnl.prop
				lnPriorRatio=nr$lnl.prior
			} else { # neither split nor merge
				newRates=adjust.rates(currRates)
				newRates.c=fix.categories(newRates)
				nRateProp=nRateProp+1
			}
			
			modCurr=lnl.BM(currRates, ape.tre, ouch.dat)
			modNew=lnl.BM(newRates, ape.tre, ouch.dat)
		}
		
		lnLikelihoodRatio = modNew$lnL - modCurr$lnL
		r=assess.lnR(heat * lnLikelihoodRatio + heat * lnPriorRatio + lnProposalRatio)
		
		if (runif(1) <= r) {			# adopt proposal
			if(currModel=="OU") {
				currAlpha <- newAlpha
				currRegimes <- newRegimes
			}
			currRates <- newRates
			currRates.c <- newRates.c
			best.lnL <- modNew$lnL
		} else {						# deny proposal			
			best.lnL <- modCurr$lnL 
		}
		
		l[i]=best.lnL
		r.nrate[i]=max(currRates.c)
		r.nregimes[i]=max(currRegimes)
		r.alpha[i]=currAlpha
		
		endparms=list(currRates, currRates.c, currAlpha, currRegimes)
		cOK=tally.mcmc.parms(startparms, endparms, cOK)
	}
# End rjMCMC
	
	summarize.run(cOK, c(nRateProp, nRcatProp, nAlphaProp, nRegimesProp), currModel)	
	if(currModel=="BM") {
		return(list(r.categories=r.nrate, lnL=l))
	} else {
		return(list(r.categories=r.nrate, r.alpha=r.alpha, r.regimes=r.nregimes, lnL=l))
	}
}


## AUXILLIARY FUNCTIONS

lnl.OU <- function(rates, alpha, regimes, apetree, ouch.dat) { # from OUCH 2.7-1:::hansen
	ouch.dat$regimes=regimes
	ouch.tre=apply.BMM.to.tree(rates, apetree)
	beta=ouch:::regime.spec(ouch.tre, ouch.dat['regimes'])
	dat=do.call(c,lapply(ouch.dat['trait'],function(y)y[ouch.tre@term]))
	res=ouch:::ou.lik.fn(ouch.tre, as.matrix(alpha), as.matrix(1), beta, dat)
	list(theta=res$coef, lnL=-0.5*res$deviance)
}

lnl.BM <- function(rates, apetree, ouch.dat) { # from OUCH 2.7-1:::brown
	ouch.tre=apply.BMM.to.tree(rates, apetree)
	dat=do.call(c,lapply(ouch.dat['trait'],function(y)y[ouch.tre@term]))
	nterm <- ouch.tre@nterm
	nchar <- 1
	w <- matrix(data=1,nrow=nterm,ncol=1)
	b <- ouch.tre@branch.times
	sols <- ouch:::glssoln(w,dat,b)
	e <- sols$residuals # residuals
	q <- t(e)%*%solve(b,e)
	v <- q/nterm
	dev <- nchar*nterm*(1+log(2*pi))+nchar*log(det(b))+nterm*log(det(v))
	list(lnL = -0.5*dev)	
}

assess.lnR <- function(lnR) {
	if(lnR > 0) {
		r=1 
	} else if(lnR < -10) {
		r=0 
	} else {
		r=exp(lnR)
	}
	r
}

decide.s.or.m <- function(categories) {
	if (max(categories) == 1) {
		return("split")
	}
	else if (max(categories) == length(categories)) {
		return("merge")
	}
	else if (runif(1) < 0.5) {
		return("split")
	}
	else return("merge")
}
	
split.or.merge <- function(values, categories, task=c("split","merge")) 
{	
	vv=values
	cc=categories
	new.vv=vv
	new.cc=cc
	nn=length(values)
	
	if(task=="merge") {
		ncat <- max(cc)
		m.cc <- sample(1:ncat)[1:2]
		new.cc[new.cc == m.cc[1]] <- m.cc[2]
		new.cc <- fix.categories(new.cc)
		or1 <- vv[cc == m.cc[1]][1]
		or2 <- vv[cc == m.cc[2]][1]
		n1 <- sum(cc == m.cc[1])
		n2 <- sum(cc == m.cc[2])
		nr <- (n1 * or1 + n2 * or2)/(n1 + n2)
		new.vv[cc == m.cc[1] | cc == m.cc[2]] <- nr
		numSplittableBeforeMerge <- sum(table(cc) > 1)
		numSplittableAfterMerge <- sum(table(new.cc) > 1)
		if (max(cc) == nn) {
			probMerge = 1
		} else {
			probMerge = 0.5
		}
		
		if (max(new.cc) == 1) {
			probSplit = 1
		} else {
			probSplit = 0.5
		}
		factor = probSplit/probMerge
		lnProposalRatio = log(factor) + log(nn) + log(ncat) - log(numSplittableAfterMerge) - log(2^(n1 + n2) - 2) - log(nr/mean(vv)) - log(n1 + n2)
		lnPriorRatio = log(Stirling2(nn, ncat)) - log(Stirling2(nn, ncat - 1))
	}
	if(task=="split") {
		ncat <- max(cc)
		while (1) {
			rcat <- round(runif(1) * ncat + 0.5)
			nc <- sum(cc == rcat)
			if (nc > 1) 
			break
		}
		while (1) {
			new <- round(runif(nc) * 2 + 0.5)
			if (length(table(new)) == 2) 
			break
		}
		new.cc[cc == rcat][new == 2] <- ncat + 1
		or <- vv[cc == rcat][1]
		n1 <- sum(new == 1)
		n2 <- sum(new == 2)
		u <- runif(1, min = -0.5 * n1 * or, max = 0.5 * n2 * or)
		nr1 <- or + u/n1
		nr2 <- or - u/n2
		new.vv[new.cc == rcat] <- nr1
		new.vv[new.cc == ncat + 1] <- nr2
		new.cc <- fix.categories(new.cc)
		numSplittableBeforeSplit <- sum(table(cc) > 1)
		numSplittableAfterSplit <- sum(table(new.cc) > 1)
		if (max(cc) == 1) {
			probSplit = 1
		} else {
			probSplit = 0.5
		}
		
		if (max(new.cc) == nn) {
			probMerge = 1
		} else {
			probMerge = 0.5
		}
		
		factor = probMerge/probSplit
		lnProposalRatio = log(factor) + log(numSplittableBeforeSplit) + log(2^(n1 + n2) - 2) + log(or/mean(vv)) + log(n1 + n2) - log(nn) - log(ncat + 1)
		lnPriorRatio = log(Stirling2(nn, ncat)) - log(Stirling2(nn, ncat + 1))
	}
	list(new.values=new.vv, new.categories=new.cc, lnl.prop=lnProposalRatio, lnl.prior=lnPriorRatio)
}

fix.categories <- function(x)
{
	f<-factor(x)
	n<-nlevels(f)
	o<-numeric(n)
	for(i in 1:n) o[i]<-which(f==levels(f)[i])[1]
	nf<-ordered(f, levels(f)[order(o)])	
	nx<-as.numeric(nf)
	nx
}

apply.BMM.to.tree <- function(rates, phy) { # SLOW STEP: convert -> deconvert tree
	phy$edge.length=phy$edge.length*rates
	ouch.tre=ape2ouch(phy)
	ouch.tre
}

adjust.value <- function(value) {
	rr=value
	rch <- runif(1, min = -0.5, max = 0.5)
	rr <- rr + rch
	rr
}

adjust.rates <- function(values) {
	vv=values
	ncat <- max(fix.categories(vv)->cc)
	rcat <- round(runif(1) * (ncat) + 0.5)
	r <- log(vv[cc == rcat][1])
	rch <- runif(1, min = -2, max = 2)
	nr <- r + rch
	vv[cc == rcat] <- exp(nr)
	vv
}

prepare.ouchdata <- function(phy, data) {
	td <- treedata(phy, data, sort = T)	
	ouch.tre <- ape2ouch(td$phy->ape.tre)		
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
	return(list(ouch.tre=ouch.tre, ouch.dat=ouch.dat, ape.tre=ape.tre, orig.dat=data))
}
			   
summarize.run<-function(cOK, cProp, currModel) {
	names(cProp)<-names(cOK)<-c("rates", "rcats", "alpha", "sRegimes")
	names(cProp)=paste("prop.",names(cProp),sep="")
	names(cOK)=paste("okay.",names(cOK),sep="")
	cat("\n")
	if(currModel=="BM") {
		print(cProp[1:2])
		print(cOK[1:2])
	} else {
		print(cProp)
		print(cOK)
	}
}

tally.mcmc.parms<-function(startparms, endparms, cOK) {
	index=array(dim=length(startparms))
	for(i in 1:length(startparms)) {
		index[i]=all(startparms[[i]]==endparms[[i]])
	}
	if(!all(index)) cOK[max(which(index==FALSE))]=cOK[max(which(index==FALSE))]+1
	cOK
}
	
Stirling2 <- function(n,m)
{
    ## Purpose:  Stirling Numbers of the 2-nd kind
    ## 		S^{(m)}_n = number of ways of partitioning a set of
    ##                      $n$ elements into $m$ non-empty subsets
    ## Author: Martin Maechler, Date:  May 28 1992, 23:42
    ## ----------------------------------------------------------------
    ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
    ## Closed Form : p.824 "C."
    ## ----------------------------------------------------------------

    if (0 > m || m > n) stop("'m' must be in 0..n !")
    k <- 0:m
    sig <- rep(c(1,-1)*(-1)^m, length= m+1)# 1 for m=0; -1 1 (m=1)
    ## The following gives rounding errors for (25,5) :
    ## r <- sum( sig * k^n /(gamma(k+1)*gamma(m+1-k)) )
    ga <- gamma(k+1)
    round(sum( sig * k^n /(ga * rev(ga))))
}

