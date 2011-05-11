## Written by LJ HARMON & AH: 2008
## Modified by JM EASTMAN: 09.2010

rjMCMC.trait<-function (phy, data, ngen=1000, sampleFreq=100, model=c("OU","BM"), probAlpha=0.05, probRegimes=0.10, probMergeSplit=0.05, probRoot=0.1, heat=1, fileBase="result") 
{
	if(length(model)>1)stop(paste("Please specify either ", sQuote("OU"), " or ", sQuote("BM")))
	
	require(ouch)
	require(geiger)
	
### prepare objects for rjMCMC
	dataList <- prepare.ouchdata(phy, data)					# check tree and data; convert into S4
	ape.tre <- currTree <- dataList$ape.tre								# S3 phylogeny
	orig.dat <- dataList$orig.dat							# S3 data
	ouch.tre <- dataList$ouch.tre							# S4 OUCH phylogeny
	ouch.dat <- dataList$ouch.dat							# S4 OUCH data
	nn <- length(ape.tre$edge.length)						# max number of rates (one per branch)
	currRates.c <- numeric(nn)								# initialize array for recording which branch belongs to which rate category
    currRates.c[] <- 1										# assign one rate to all branches initially
    currRates <- numeric(nn)								# initialize array for rate estimates branchwise
    currRates[] <- 1										# initialize BMMrate to be small
	currAlpha <- 0.001										# initialize alpha
	currRegimes <- numeric(nrow(ouch.dat))					# initial set of selective regimes
	currRegimes[] <- 1										# initialize regime to be global
	currTheta <- 0											# initialize trait optima under OU
	currRoot <- mean(ouch.dat['trait'], na.rm=TRUE)
	currModel = model										# initialize family of models explored
    nRateProp = 0											# proposal count: update rate parameter(s)
	nRcatProp = 0											# proposal count: update rate categories
	nRootProp = 0
	nAlphaProp = 0											# proposal count: update alpha parameter
	nRegimesProp = 0										# proposal count: update selective regimes
	nThetaProp = 0
	nRateOK = 0												# proposal count: successful rate changes
	nRcatOK = 0												# proposal count: successful rate categories updates
	nRootOK = 0
	nAlphaOK = 0											# proposal count: successful alpha updates
	nRegimesOK = 0											# proposal count: successful regime updates
	nThetaOK = 0
	l <- array(dim=c(ngen,4))								# array for storing output: size is number of MCMC generations
    l[] <- NA												# initialize array for storing lnLs
	
	cOK=c(													# list of proposal tallies
		  nRateOK, 
		  nRcatOK, 
		  nRootOK, 
		  nAlphaOK, 
		  nRegimesOK, 
		  nThetaOK
		  )		
	
	errorLog=file(paste(getwd(),paste(currModel, "rjTraitErrors.log",sep="."),sep="/"),open='w+')
	generate.error.message(i=NULL, modCurr=NULL, modNew=NULL, lnR=NULL, errorLog=errorLog, init=TRUE)
	runLog=file(paste(getwd(),paste(currModel, "rjTrait.log",sep="."),sep="/"),open='w+')
	generate.log(bundled.parms=NULL, currModel, file=runLog, init=TRUE)
	
### Begin rjMCMC
    for (i in 1:ngen) {
		error.flag=FALSE
        lnLikelihoodRatio <- lnProposalRatio <- lnPriorRatio <- 0
		startparms=list(rates=currRates, cats=currRates.c, root=currRoot, alph=currAlpha, reg=currRegimes, theta=currTheta)
		
## OU IMPLEMENTATION ## 
		if(currModel=="OU") {
			if (runif(1) < (2 * probMergeSplit)) {
				if(runif(1)<0.5) { # adjust regimes
					nr=split.or.merge(currTheta, currRegimes, decide.s.or.m(currRegimes, regimes=TRUE)->s.m, regimes=TRUE)
					newRegimes=nr$new.categories
					newTheta=nr$new.values
					newRates=currRates
					newRates.c=currRates.c
					nRegimesProp=nRegimesProp+1
				} else { # adjust rate categories
					nr=split.or.merge(currRates, currRates.c, decide.s.or.m(currRates.c)->s.m)
					newRates=nr$new.values
					newRates.c=nr$new.categories
					newRegimes=currRegimes
					newTheta=currTheta
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
					newRegimes=currRegimes
					newTheta=currTheta
					nAlphaProp=nAlphaProp+1
				} else { 
					if(runif(1)<0.5) {# adjust rates
						newRates=adjust.rates(currRates)
						newRates.c=fix.categories(newRates)
						newAlpha=currAlpha
						newRegimes=currRegimes
						newTheta=currTheta
						nRateProp=nRateProp+1
					} else { # adjust theta
						newTheta=sapply(currTheta, adjust.value)
						newAlpha=currAlpha
						newRates=currRates
						newRates.c=currRates.c
						newRegimes=currRegimes
						nThetaProp=nThetaProp+1
					}
				}
			}
			
			currTree=apply.BMM.to.apetree(ape.tre, currRates)
			newTree=apply.BMM.to.apetree(ape.tre, newRates)
			modCurr=new.ou.lik.fn(currRates, currAlpha, currRegimes, currTheta, currTree, ouch.dat)
			modNew=try(new.ou.lik.fn(newRates, newAlpha, newRegimes, newTheta, newTree, ouch.dat),silent=TRUE)
		} else if(currModel=="BM"){
## BM IMPLEMENTATION ##
			if (runif(1) < (2 * probMergeSplit)) {
				nr=split.or.merge(currRates, currRates.c, decide.s.or.m(currRates.c)->s.m)
				newRates=nr$new.values
				newRates.c=nr$new.categories
				newRoot=currRoot
				nRcatProp=nRcatProp+1
				lnProposalRatio=nr$lnl.prop
				lnPriorRatio=nr$lnl.prior
			} else { 
				if(runif(1)<probRoot) { # adjust root
					newRoot=adjust.value(currRoot)
					newRates=currRates
					newRates.c=currRates.c
					nRootProp=nRootProp+1
				} else { # adjust rates
					newRates=adjust.rates(currRates)
					newRates.c=fix.categories(newRates)
					newRoot=currRoot
					nRateProp=nRateProp+1
				}
			}
			
			currTree=apply.BMM.to.apetree(ape.tre, currRates)
			newTree=apply.BMM.to.apetree(ape.tre, newRates)
			modCurr=new.bm.lik.fn(currRates, currRoot, currTree, ouch.dat)
			modNew=try(new.bm.lik.fn(newRates, newRoot, newTree, ouch.dat),silent=TRUE)
		}
	
		if(!inherits(modNew, "try-error")) {
			lnLikelihoodRatio = modNew$lnL - modCurr$lnL
		} else {
			lnLikelihoodRatio = -Inf
			modNew$lnL=NaN
		}
		r=assess.lnR((heat * lnLikelihoodRatio + heat * lnPriorRatio + lnProposalRatio)->lnR)
		
		if (runif(1) <= r$r) {			# adopt proposal
			if(currModel=="OU") {
				currAlpha <- newAlpha
				currRegimes <- newRegimes
				currTheta <- newTheta
			} else {
				currRoot <- newRoot
			}
			currRates <- newRates
			currRates.c <- newRates.c
			curr.lnL <- modNew$lnL
			currTree <- newTree
		} else {						# deny proposal			
			curr.lnL <- modCurr$lnL 
		}
		if(r$error) generate.error.message(i, modCurr, modNew, lnR, errorLog)
		
		l[i,]<-c(curr.lnL, max(currRates.c), currAlpha, max(currRegimes))
		endparms=list(rates=currRates, cats=currRates.c, root=currRoot, alph=currAlpha, reg=currRegimes, theta=currTheta)
		cOK=tally.mcmc.parms(startparms, endparms, cOK)
		
		bundled.parms=list(gen=i, mrate=mean(currRates), cats=max(currRates.c), root=currRoot, alpha=currAlpha, reg=max(currRegimes), theta=mean(currTheta), lnL=curr.lnL)
		if(i%%sampleFreq==0) generate.log(bundled.parms, currModel, file=runLog)
		print(currRates)
	}
	
# End rjMCMC
	close(errorLog)
	close(runLog)
	summarize.run(cOK, (cProp<-c(nRateProp, nRcatProp, nRootProp, nAlphaProp, nRegimesProp, nThetaProp)), currModel)	
	if(currModel=="BM") {
		return(list(r.cat=l[,2], lnL=l[,1]))
	} else {
		return(list(r.cat=l[,2], r.alpha=l[,3], r.regimes=l[,4], lnL=l[,1]))
	}
}


## AUXILIARY FUNCTIONS


new.ou.lik.fn <- function(rates, alpha, regimes, theta, ape.tre, ouch.dat) { # from OUCH 2.7-1:::hansen()
	ouch.tre=ape2ouch(ape.tre, scale=TRUE)
	ouch.dat$regimes=as.factor(regimes)
	alpha=as.matrix(sqrt(alpha))
	beta=ouch:::regime.spec(ouch.tre, ouch.dat['regimes'])
	dat=do.call(c,lapply(ouch.dat['trait'],function(y)y[ouch.tre@term]))
	n <- length(dat)
	ev <- eigen(alpha,symmetric=TRUE)
	w <- .Call(ouch:::ouch_weights,object=ouch.tre,lambda=ev$values,S=ev$vectors,beta=beta)
	v <- .Call(ouch:::ouch_covar,object=ouch.tre,lambda=ev$values,S=ev$vectors,sigma.sq=1)
	y <- matrix(theta,ncol=1)
	e <- w%*%y-dat
	dim(y) <- dim(y)[1]
	dim(e) <- n
	q <- e%*%solve(v,e)
	det.v <- determinant(v,logarithm=TRUE)
	if (det.v$sign!=1)stop("ou.lik.fn error: non-positive determinant",call.=FALSE)
	dev <- n*log(2*pi)+as.numeric(det.v$modulus)+q[1,1]
	list(
		 theta=theta,
		 weight=w,
		 vcov=v,
		 resids=e,
		 lnL=-0.5*dev
		 )
}

new.bm.lik.fn <- function(rates, root, ape.tre, ouch.dat) { # from OUCH 2.7-1:::brown()
	ouch.tre=ape2ouch(ape.tre, scale=TRUE)
	dat=do.call(c,lapply(ouch.dat['trait'],function(y)y[ouch.tre@term]))
	nterm <- ouch.tre@nterm
	nchar <- 1
	w <- matrix(data=1,nrow=nterm,ncol=1)
	b <- ouch.tre@branch.times
	y <- matrix(root)
	e <- w%*%y-dat
	q <- t(e)%*%solve(b,e)
	v <- q/nterm
	dev <- nchar*nterm*(1+log(2*pi))+nchar*log(det(b))+nterm*log(det(v))
	list(
		 root=root,
		 lnL = -0.5*dev
		 )	
}

assess.lnR <- function(lnR) {
	if(is.na(lnR) || abs(lnR)==Inf) {
		error=TRUE
		r=0
	} else {
		if(lnR < -100) {
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

decide.s.or.m <- function(categories, regimes=FALSE) {
	if (max(categories) == 1) {
		return("split")
	}
	else if (max(categories) == length(categories)) {
		return("merge")
	}
#	else if (regimes && (((length(categories)+1)/2)==max(categories))) {
# limits number of selective regimes to number of tips
#		return("merge")
#	}
	else if (runif(1) < 0.5) {
		return("split")
	}
	else return("merge")
}
	
split.or.merge <- function(values, categories, task=c("split","merge"), regimes=FALSE) 
{	
# updated to adjust for an additional theta when regimes are split (or one fewer theta if regimes are merged)
# implement Drummon and Suchard 2010 method, adapting local relaxed clocks for local theta and BMMrates - draw changes in BMMrates from truncated Poisson
	cc=categories

	if(regimes) vv=values[cc] else vv=values
	new.vv=vv
	new.cc=cc
	nn=length(vv)
	
	if(task=="merge") {
		ncat <- max(cc)
		m.cc <- sample(1:ncat)[1:2] # grab two categories
		new.cc[new.cc == m.cc[1]] <- m.cc[2]	# merge selected categories
		new.cc <- fix.categories(new.cc)
		or1 <- vv[cc == m.cc[1]][1] # value affected
		or2 <- vv[cc == m.cc[2]][1]	# value affected
		n1 <- sum(cc == m.cc[1]) # no. entries affected
		n2 <- sum(cc == m.cc[2]) # no. entries affected
		nr <- (n1 * or1 + n2 * or2)/(n1 + n2) # find new value for merged category
		new.vv[cc == m.cc[1] | cc == m.cc[2]] <- nr # new value assigned to merged classes
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
# AS WRITTEN by LJ HARMON
		lnProposalRatio = log(factor) + log(nn) + log(ncat) - log(numSplittableAfterMerge) - log(2^(n1 + n2) - 2) - log(nr/mean(vv)) - log(n1 + n2)
		lnPriorRatio = log(Stirling2(nn, ncat)) - log(Stirling2(nn, ncat - 1))
	}
	if(task=="split") {
		ncat <- max(cc)
		while (1) { # find non-unique category
			rcat <- round(runif(1) * ncat + 0.5) # class to split
			nc <- sum(cc == rcat) 
			if (nc > 1) 
			break
		}
		while (1) { # 
			new <- round(runif(nc) * 2 + 0.5)
			if (length(table(new)) == 2) 
			break
		}
		new.cc[cc == rcat][new == 2] <- ncat + 1
		or <- vv[cc == rcat][1] # old rate
		n1 <- sum(new == 1)	# no. elements affected
		n2 <- sum(new == 2)	# no. elements affected
		lims <- c(-0.5 * n1 * or, 0.5 * n2 * or)
		lims <- lims[order(lims)]
		u <- runif(1, min = lims[1], max = lims[2])
		nr1 <- or + u/n1 # new rate 1
		nr2 <- or - u/n2 # new rate 2
		new.vv[new.cc == rcat] <- nr1 # assign new rate 1
		new.vv[new.cc == ncat + 1] <- nr2 # assign new rate 2
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
		
# AS WRITTEN by LJ HARMON
		lnProposalRatio = log(factor) + log(numSplittableBeforeSplit) + log(2^(n1 + n2) - 2) + log(or/mean(vv)) + log(n1 + n2) - log(nn) - log(ncat + 1)
		lnPriorRatio = log(Stirling2(nn, ncat)) - log(Stirling2(nn, ncat + 1))
	}
	if(regimes) new.vv=unique(new.vv[order(new.cc)])
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

assign.edgelengths<-function(ouchtree) {
	ouchtree=as(ouchtree, "data.frame")
	ouchtree$edges=NA
	for(i in 2:nrow(ouchtree)) {
		ouchtree$edges[i]=ouchtree$times[i]-(ouchtree[match(ouchtree[i,'ancestors'], ouchtree[,'nodes']),'times'])
	}
	ouchtree$edges[1]=0
	ouchtree
}

apply.BMM.to.ouchtree<-function(ouchtree, rates) {
## FIXME: ensure that the tree maintains structure before and after this function is called
## VERY SLOW 
## NOT CURRENTLY implemented
	ot=as(ouchtree, "data.frame")
	ot=assign.edgelengths(ot)
	ot$edges.new=c(0, rates)*ot$edges
	ot=ot[order(ot$times),]
	ot$times.new=NA
	ot$times.new[1]=0
	for(i in 2:nrow(ot)) {
		ot$times.new[i]=ot[match(ot[i,'ancestors'], ot[,'nodes']),'times.new']+ot[i,'edges.new']
	}
	ot$times=ot$times.new
	return(with(ot, ouchtree(nodes=nodes, ancestors=ancestors, times=times, labels=labels)))
}

apply.BMM.to.apetree<-function(apetree, rates) {
## FASTER than apply.BMM.to.ouchtree
	apetree$edge.length=apetree$edge.length*rates
	return(apetree)
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
	ouch.tre <- ape2ouch(td$phy->ape.tre, scale=TRUE)		
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
	names(cProp)<-names(cOK)<-c("rates", "rcats", "root", "alpha", "sRegimes", "theta")
	names(cProp)=paste("prop.",names(cProp),sep="")
	names(cOK)=paste("okay.",names(cOK),sep="")
	cat("\n")
	if(currModel=="BM") {
		print(cProp[1:3])
		print(cOK[1:3])
	} else {
		print(cProp[c(1,2,4,5,6)])
		print(cOK[c(1,2,4,5,6)])
	}
}

generate.log<-function(bundled.parms, currModel, file, init=FALSE) {
#	bundled.parms=list(gen=i, mrate=mean(currRates), cats=max(currRates.c), root=currRoot, alpha=currAlpha, reg=max(currRegimes), theta=mean(currTheta), lnL=curr.lnL)
	res=bundled.parms

	if(init) {
		if(currModel=="BM") msg=paste("gen", "model", "mean.rate", "rates", "root", "lnL", sep="\t") else msg=paste("gen", "model", "mean.rate", "rates", "alpha", "regimes", "mean.theta", "lnL", sep="\t")
	} else {
		if(currModel=="BM") {
			msg<-paste(res$gen, sQuote(currModel), sprintf("%.3f", res$mrate), res$cats, sprintf("%.3f", res$root), sprintf("%.3f", res$lnL),sep="\t")
		} else {
			msg<-paste(res$gen, sQuote(currModel), sprintf("%.3f", res$mrate), res$cats, sprintf("%.4f", res$alpha), res$reg, sprintf("%.3f", res$theta), sprintf("%.3f", res$lnL),sep="\t")
		}
	}
	write(msg, file=file, append=TRUE)
}

generate.error.message<-function(i, modCurr, modNew, lnR, errorLog, init=FALSE) {
	if(init) {
		initmsg = write(paste("gen", "curr.lnL", "new.lnL", "lnR", sep="\t"), file=errorLog)
	} else {
		write(paste(i, sprintf("%.3f", modCurr$lnL), sprintf("%.3f", modNew$lnL), sprintf("%.3f", lnR), sep="\t"), file=errorLog)
	}
}

tally.mcmc.parms<-function(startparms, endparms, cOK) {
	index=array(dim=(length(startparms)-1))
	for(i in 1:(length(startparms)-1)) {
		index[i]=all(startparms[[i]]==endparms[[i]])
	}
	if(!all(index)) {
		cOK[max(which(index==FALSE))]=cOK[max(which(index==FALSE))]+1
	}
	if(length(startparms$theta) == length(endparms$theta)) {
		if(!all(sort(startparms$theta)==sort(endparms$theta))) {
			cOK[6]=cOK[6]+1
		}
	}
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

