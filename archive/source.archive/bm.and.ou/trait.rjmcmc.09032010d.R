## Written by LJ HARMON, AH: 2008
## Updated by JM EASTMAN: 09.2010

rjMCMC.trait<-function (phy, data, ngen = 1000, changeModel = 0.05, probAlpha=0.05, probRegimes=0.10, probMergeSplit = 0.05, probRoot = 0.1, sampleFreq = 100, ylim = c(-200, 0), heat=1, fileBase = "result") 
{
	errorNum <- 0 
	errorLog <- paste(fileBase, 'errors.log', sep = '')
	cat(c('Error log for run beginning', date(), '\n'), file = errorLog, append = F, sep = ' ')
	cat(c('--------------------------------------------------------', '\n\n'), file = errorLog, append = T)
	cat(c('generation', '\t', 'numCat', '\t', 'rateMin', '\t', 'rateMax', '\t', 'cumulativePercentage', '\n'), file = paste(errorLog, 'stats.txt', sep = ''), append = F, sep = '')
	
	if (fileBase != "") {
        f1 <- paste(fileBase, "RateCat.txt", sep = "")
        f2 <- paste(fileBase, "Rates.txt", sep = "")
    } else {
        f1 <- ""
        f2 <- ""
    }
    cat("", file = f1)
    cat("", file = f2)
	
### construct a dataset interpretable by OUCH
	require(ouch)
	require(geiger)
	td <- treedata(phy, data, sort = T)	
	ouch.tre <- ape2ouch(td$phy->ape.tre)		
	d <- data
	ouch.dat=data.frame(as.numeric(ouch.tre@nodes), as.character(ouch.tre@nodelabels))
	names(ouch.dat)=c("nodes","labels")
	ouch.dat[ouch.dat[,2]=="",2]=NA
	ouch.dat$trait=NA
	mM=match(ouch.dat$labels,names(d))
	for(i in 1:nrow(ouch.dat)){
		if(!is.na(mM[i])){
			ouch.dat$trait[i]=d[mM[i]]	
		}
	}
	dat=do.call(c,lapply(ouch.dat['trait'],function(y)y[ouch.tre@term]))

###
	
### prepare objects for rjMCMC
    nn <- length(td$phy$edge.length)		# max number of rates (one per branch)
	currRates.c <- numeric(nn)				# initialize array for recording which branch belongs to which rate category
    currRates.c[] <- 1						# assign one rate to all branches initially
    currRates <- numeric(nn)					# initialize array for rate estimates branchwise
    currRates[] <- init <- 0.001				# initialize BMMrate to be small
    currRootMean <- mean(td$data[, 1])		# determine an ancestral state (by averaging known data)
#	currTheta <- numeric(nn)				# max number of optima (one per branch)
#	currTheta[] <- mean(td$data[,1])		# initialize optima
	currAlpha <- 0.001						# initialize alpha
	currRegimes <- numeric(nrow(ouch.dat))	# initial set of selective regimes
	currRegimes[] <- 1						# initialize regime to be global
	currModel <- "OU"
#	currModel <- sample(c("BM","OU"),1)		# initialize model of trait evolution
	nMergeProp = 0							# proposal count: merge two rate parameters
    nSplitProp = 0							# proposal count: draw new rate parameter
    nRateProp = 0							# proposal count: update rate parameter(s)
#	nThetaProp = 0							# proposal count: update theta parameter(s)
	nAlphaProp = 0							# proposal count: update alpha parameter
    nMergeOK = 0							# proposal count: successful merges
    nSplitOK = 0							# proposal count: successful splits
    nRateOK = 0								# proposal count: successful rate changes
#	nThetaProp = 0							# proposal count: successful theta updates
	nAlphaProp = 0							# proposal count: successful alpha updates
	l <- array(dim=ngen)					# number of MCMC generations
    l[] <- NA								# initialize array for storing lnLs
###
	
# Begin rjMCMC
    for (i in 1:ngen) {
        lnLikelihoodRatio <- lnProposalRatio <- lnPriorRatio <- 0
		
		if(runif(1) < changeModel) {
			if(currModel=="BM") currModel="OU" else currModel="BM"
		}
		
## OU IMPLEMENTATION ## thetas currently optimized by ouch:::glssoln(); selective regimes are choosable (see Beta)
		if(currModel=="OU") { # proposal mechanism appears working
			if (runif(1) < (2 * probMergeSplit)) {
				adj.regimes <- adj.rates <- FALSE
				if(runif(1) < 0.5) adj.regimes=TRUE else adj.rates=TRUE
				if(adj.regimes) {
					newRegimes=split.or.merge(currRegimes, currRegimes, decide.s.or.m(currRegimes))$new.categories
					newRates=currRates
					newRates.c=currRates.c
				} else { # adjust rate categories
					nr=split.or.merge(currRates, currRates.c, decide.s.or.m(currRates.c))
					newRates=nr$new.values
					newRates.c=nr$new.categories
					newRegimes=currRegimes
				}
				newAlpha=currAlpha
			} else { # neither split nor merge
				if(runif(1)<probAlpha) { # adjust alpha
					newAlpha=adjust.rates(currAlpha)
					newRates=currRates
					newRates.c=currRates.c
				} else { # adjust rates
					newRates=adjust.rates(currRates)
					newRates.c=fix.categories(newRates)
					newAlpha=currAlpha
				}
				newRegimes=currRegimes
			}			
			modCurr=lnl.OU(currRates, currAlpha, currRegimes, ape.tre, ouch.dat)
			modNew=lnl.OU(newRates, newAlpha, newRegimes, ape.tre, ouch.dat)
			
			lnLikelihoodRatio = modNew$lnL - modCurr$lnL
			r=assess.lnR(heat * lnLikelihoodRatio + heat * lnPriorRatio + lnProposalRatio)
						
# adopt proposed update; message("YES")
			if (runif(1) <= r) {
				currRates <- newRates
				currRates.c <- newRates.c
				currAlpha <- newAlpha
				currRegimes <- newRegimes
				best.lnL <- modNew$lnL
			} else {					
				best.lnL <- modCurr$lnL  
			}
		} else if(currModel=="BM"){ # BM
			nterm <- ouch.tre@nterm
			nchar <- 1
			w <- matrix(data=1,nrow=nterm,ncol=1)
			b <- ouch.tre@branch.times
			sols <- ouch:::glssoln(w,dat,b)
			theta <- sols$coeff # optima
			e <- sols$residuals # residuals
			q <- t(e)%*%solve(b,e)
			v <- q/nterm
			dev <- nchar*nterm*(1+log(2*pi))+nchar*log(det(b))+nterm*log(det(v))
			lnlCurr<-bm.log.lik<- -0.5*dev
		}
		l[i]=best.lnL
    }
# End rjMCMC
	return(list(alpha=currAlpha, rates=max(currRates.c), regimes=max(currRegimes), l))
}

lnl.OU<-function(rates, alpha, regimes, apetree, ouch.dat) {
	ouch.dat$regimes=regimes
	ouchtree = apply.rates.to.tree(rates, apetree)
	beta = ouch:::regime.spec(ouchtree, ouch.dat['regimes'])
	dat=do.call(c,lapply(ouch.dat['trait'],function(y)y[ouchtree@term]))
	res=ouch:::ou.lik.fn(ouchtree, as.matrix(alpha), as.matrix(1), beta, dat)
	list(theta=res$coef, lnL=-0.5*res$deviance)
}

assess.lnR<-function(lnR) {
	if(lnR > 0) {
		r=1 
	} else if(lnR < -100) {
		r=0 
	} else {
		r=exp(lnR)
	}
	r
}

decide.s.or.m<-function(categories) {
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
	
split.or.merge<-function(values, categories, task=c("split","merge")) 
# generalized: to work for either BMM rates or for OU optima
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
	list(new.values=new.vv, new.categories=new.cc, lnProp=lnProposalRatio, lnPrior=lnPriorRatio)
}

fix.categories<-function(x)
# assigns indicator variables to shared rate classes amongst branches
{
	f<-factor(x)
	n<-nlevels(f)
	o<-numeric(n)
	for(i in 1:n) o[i]<-which(f==levels(f)[i])[1]
	nf<-ordered(f, levels(f)[order(o)])	
	nx<-as.numeric(nf)
	nx
}

apply.rates.to.tree<-function(rates, phy) {
	phy$edge.length=phy$edge.length*rates
	ouch.tre=ape2ouch(phy)
	ouch.tre
}

adjust.value<-function(value) {
	rr=value
	rch <- runif(1, min = -0.5, max = 0.5)
	rr <- rr + rch
	rr
}

adjust.rates<-function(values) {
	vv=values
	ncat <- max(fix.categories(vv)->cc)
	rcat <- round(runif(1) * (ncat) + 0.5)
	r <- log(vv[cc == rcat][1])
	rch <- runif(1, min = -2, max = 2)
	nr <- r + rch
	vv[cc == rcat] <- exp(nr)
	vv
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

