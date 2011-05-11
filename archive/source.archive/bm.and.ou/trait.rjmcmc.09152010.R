## Written by LJ HARMON & AH: 2008
## Modified by JM EASTMAN: 09.2010

rjMCMC.trait<-function (phy, data, ngen=1000, sampleFreq=100, model=c("OU","BM"), probAlpha=0.05, probMergeSplit=0.05, probRoot=0.1, heat=1, fileBase="result") 
{
	if(length(model)>1)stop(paste("Please specify either ", sQuote("OU"), " or ", sQuote("BM")))
## MORE ERROR CHECKING NEEDED HERE ## 
	
	require(ouch)
	require(geiger)
	
### prepare objects for rjMCMC
	cur.model		<- model								# initialize family of models explored

	dataList		<- prepare.ouchdata(phy, data)			# check tree and data; convert into S4

	ape.tre			<- cur.tre <- dataList$ape.tre			# S3 phylogeny
	orig.dat		<- dataList$orig.dat					# S3 data	
	ouch.tre		<- dataList$ouch.tre					# S4 OUCH phylogeny
	ouch.dat		<- dataList$ouch.dat					# S4 OUCH data
	
	nn				<- length(ape.tre$edge.length)			# max number of rates (one per branch)
	init.rate		<- fitContinuous(ape.tre,orig.dat)
    cur.rates		<- rep(init.rate[[1]]$beta,nn)			# initialize array for rate estimates branchwise
	cur.delta.rates	<- rep(0,nn)
	cur.root		<- mean(orig.dat)
	cur.alpha		<- 0.001								# initialize alpha
	cur.theta		<- rep(0,nn)							# initialize trait optima under OU
	cur.delta.theta	<- rep(0,nn)							# initialize local 'rates'
	cur.regimes		<- as.factor(c(sample(cur.theta,1),cur.theta))				# initial set of selective regimes
	cur.root.theta	<- sample(cur.theta,1)

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
	
	errorLog=file(paste(getwd(),paste(cur.model, "rjTraitErrors.log",sep="."),sep="/"),open='w+')
	generate.error.message(i=NULL, mod.cur=NULL, mod.new=NULL, lnR=NULL, errorLog=errorLog, init=TRUE)
	runLog=file(paste(getwd(),paste(cur.model, "rjTrait.log",sep="."),sep="/"),open='w+')
	generate.log(bundled.parms=NULL, cur.model, file=runLog, init=TRUE)
	
### Begin rjMCMC
    for (i in 1:ngen) {
        lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
 startparms=list(rates=cur.rates, cats=sum(cur.delta.rates)+1, root=cur.root, alph=cur.alpha, reg=length(table(cur.regimes)), theta=mean(cur.theta))
		
## OU IMPLEMENTATION ## 
		if(cur.model=="OU") {
			if (runif(1) < (2 * probMergeSplit)) {
				if(runif(1)<0.5) {				# adjust theta categories
					nr=split.or.merge(cur.delta.theta, cur.theta, ape.tre, theta=TRUE)
					new.alpha=cur.alpha
					new.root.theta=cur.root.theta
					new.theta=nr$new.values ##
					new.delta.theta=nr$new.delta ##
					nrt=assign.root.theta(new.theta, ape.tre) ## 
					new.regimes=as.factor(c(nrt,new.theta)) ##
					new.rates=cur.rates 
					new.delta.rates=cur.delta.rates 
					nRegimesProp=nRegimesProp+1
					print("try.theta")
				} else {						# adjust rate categories
					nr=split.or.merge(cur.delta.rates, cur.rates, ape.tre, theta=FALSE)
					new.alpha=cur.alpha
					new.root.theta=cur.root.theta
					new.theta=cur.theta
					new.delta.theta=cur.delta.theta
					new.regimes=cur.regimes
					new.rates=nr$new.values ##
					new.delta.rates=nr$new.delta ##
					nRcatProp=nRcatProp+1
					print("try.rates")
				}
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
			} else { # neither split nor merge
				if(runif(1)<probAlpha) {		# adjust alpha
					new.alpha=adjust.rate(cur.alpha) ##
					new.root.theta=cur.root.theta
					new.theta=cur.theta
					new.delta.theta=cur.delta.theta
					new.regimes=cur.regimes 
					new.rates=cur.rates
					new.delta.rates=cur.delta.rates
					nAlphaProp=nAlphaProp+1
				} else { 
					if(runif(1)<0.5) {			# adjust rates
						new.alpha=cur.alpha
						new.root.theta=cur.root.theta
						new.theta=cur.theta
						new.delta.theta=cur.delta.theta
						new.regimes=cur.regimes 
						new.rates=adjust.rate(cur.rates) ##
						new.delta.rates=cur.delta.rates
						nRateProp=nRateProp+1
					} else {					# adjust theta
						new.alpha=cur.alpha 
						new.root.theta=cur.root.theta
						new.theta=adjust.value(cur.theta) ##
						new.delta.theta=cur.delta.theta
						nrt=assign.root.theta(new.theta, ape.tre) ## 
						new.regimes=as.factor(c(nrt,new.theta)) ##
						new.rates=cur.rates
						new.delta.rates=cur.delta.rates
						nThetaProp=nThetaProp+1				
					}
				}
			}
			
			cur.tre=apply.BMM.to.apetree(ape.tre, cur.rates)
			new.tre=apply.BMM.to.apetree(ape.tre, new.rates)
			mod.cur=new.ou.lik.fn(cur.rates, cur.alpha, cur.regimes, cur.theta, cur.tre, ouch.dat)
			mod.new=new.ou.lik.fn(new.rates, new.alpha, new.regimes, new.theta, new.tre, ouch.dat)
#			mod.new=try(new.ou.lik.fn(new.rates, new.alpha, new.regimes, new.theta, new.tre, ouch.dat),silent=TRUE)
		} else if(cur.model=="BM"){
## BM IMPLEMENTATION ##
			if (runif(1) < (2 * probMergeSplit)) {
				nr=split.or.merge(cur.delta.rates, cur.rates, ape.tre, theta=FALSE)
				new.rates=nr$new.values
				new.delta.rates=nr$new.delta
				new.root=cur.root
				nRcatProp=nRcatProp+1
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
			} else { 
				if(runif(1)<probRoot) { # adjust root
					new.root=adjust.value(cur.root)
					new.rates=cur.rates
					new.delta.rates=cur.delta.rates
					nRootProp=nRootProp+1
				} else { # adjust rates
					new.rates=adjust.rate(cur.rates)
					new.delta.rates=cur.delta.rates
					new.root=cur.root
					nRateProp=nRateProp+1
				}
			}
			
			cur.tre=apply.BMM.to.apetree(ape.tre, cur.rates)
			new.tre=apply.BMM.to.apetree(ape.tre, new.rates)
			mod.cur=new.bm.lik.fn(cur.rates, cur.root, cur.tre, ouch.dat)
			mod.new=new.bm.lik.fn(new.rates, new.root, new.tre, ouch.dat)
#			mod.new=try(new.bm.lik.fn(new.rates, new.root, new.tre, ouch.dat),silent=TRUE)
		}
	
		if(!inherits(mod.new, "try-error")) {
			lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
		} else {
			lnLikelihoodRatio = -Inf
			mod.new$lnL=NaN
		}
		r=assess.lnR((heat * lnLikelihoodRatio + heat * lnPriorRatio + lnHastingsRatio)->lnR)

		if (runif(1) <= r$r) {			# adopt proposal
			if(cur.model=="OU") {
				decision="adopt"
				cur.alpha <- new.alpha
				cur.root.theta <- new.root.theta
				cur.regimes <- new.regimes
				cur.theta <- new.theta
				cur.delta.theta <- new.delta.theta
			} else {
				decision="adopt"
				cur.root <- new.root
			}
			cur.rates <- new.rates
			cur.delta.rates <- new.delta.rates
			curr.lnL <- mod.new$lnL
			cur.tre <- new.tre
		} else {						# deny proposal	
			decision="reject"
			curr.lnL <- mod.cur$lnL 
		}
		
		if(i==1) cat(paste("gen", "cur.lnL", "new.lnL", "lnR", sQuote(decision), sep="\t\t"),"\n")
		cat(paste(i, sprintf("%.2f", mod.cur$lnL), sprintf("%.2f", mod.new$lnL), sprintf("%.2f", lnR), sQuote(decision), sep="\t\t"),"\n")

		if(r$error) generate.error.message(i, mod.cur, mod.new, lnR, errorLog)
		
		l[i,]<-c(curr.lnL, sum(cur.delta.rates)+1, cur.alpha, length(table(cur.regimes)))
		endparms=list(rates=cur.rates, cats=sum(cur.delta.rates)+1, root=cur.root, alph=cur.alpha, reg=length(table(cur.regimes)), theta=mean(cur.theta))
		cOK=tally.mcmc.parms(startparms, endparms, cOK)
		
		bundled.parms=list(gen=i, mrate=exp(mean(log(cur.rates))), cats=sum(cur.delta.rates)+1, root=cur.root, alpha=cur.alpha, reg=length(table(cur.regimes)), theta=mean(cur.theta), lnL=curr.lnL)
		if(i%%sampleFreq==0) generate.log(bundled.parms, cur.model, file=runLog)
	}
	
# End rjMCMC
	close(errorLog)
	close(runLog)
	summarize.run(cOK, (cProp<-c(nRateProp, nRcatProp, nRootProp, nAlphaProp, nRegimesProp, nThetaProp)), cur.model)	
	if(cur.model=="BM") {
		return(list(r.cat=l[,2], lnL=l[,1]))
	} else {
		return(list(r.cat=l[,2], r.alpha=l[,3], r.regimes=l[,4], lnL=l[,1]))
	}
}


## AUXILIARY FUNCTIONS


new.ou.lik.fn <- function(rates, alpha, regimes, theta, ape.tre, ouch.dat) { # from OUCH 2.7-1:::hansen()
# mod 09.13.2010 JM Eastman
	ouch.tre=ape2ouch(ape.tre, scale=FALSE)
	ouch.dat$regimes=as.factor(regimes)
	theta.tmp=as.numeric(levels(as.factor(theta)))
	theta=theta.tmp[order(match(theta.tmp,theta))]
	alpha.ouch=as.matrix(sqrt(alpha))
	beta=ouch:::regime.spec(ouch.tre, ouch.dat['regimes'])
	dat=do.call(c,lapply(ouch.dat['trait'],function(y)y[ouch.tre@term]))
	names(dat)=ouch.dat$labels[Ntip(ape.tre):nrow(ouch.dat)]
	n <- length(dat)
	otree.df=as(ouch.tre,"data.frame")
	ordering=otree.df[n:nrow(otree.df),'labels']
	ev <- eigen(alpha.ouch,symmetric=TRUE)
	w <- .Call(ouch:::ouch_weights,object=ouch.tre,lambda=ev$values,S=ev$vectors,beta=beta)
	v <- geiger:::ouMatrix(vcv.phylo(ape.tre), alpha)
	v.match <- match(row.names(v),ordering)
	v <- v[v.match,v.match]
		## GEIGER call does not assume ultrametricity (and thus allows different 'rates' under BM)
		# -old-	v <- .Call(ouch:::ouch_covar,object=ouch.tre,lambda=ev$values,S=ev$vectors,sigma.sq=1)
	
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
# mod 09.13.2010 JM Eastman
	tree=ape2ouch(ape.tre, scale=FALSE)
	data=ouch.dat['trait']
	
	if (!is(tree, "ouchtree")) {stop(sQuote("tree"), " must be an object of class ", sQuote("ouchtree"))}
    if (is.data.frame(data)) {
        nm <- rownames(data)
        data <- lapply(as.list(data), function(x) {names(x) <- nm; x})
    }
    if (is.numeric(data)) {
        nm <- deparse(substitute(data))[1]
        data <- list(data)
        names(data) <- nm
    }
    if (is.list(data)) {
        if (any(sapply(data, class) != "numeric") || any(sapply(data, length) != tree@nnodes)) {stop(sQuote("data"), " vector(s) must be numeric, with one entry per node of the tree")}
        if (any(sapply(data, function(x) (is.null(names(x))) || (!setequal(names(x), tree@nodes))))) {stop(sQuote("data"), " vector names (or data-frame row names) must match node names of ", sQuote("tree"))}
        for (xx in data) {
            no.dats <- which(is.na(xx[tree@nodes[tree@term]]))
            if (length(no.dats) > 0) {stop("missing data on terminal node(s): ", paste(tree@nodes[tree@term[no.dats]], collapse = ","))}
        }
    }
    nterm <- tree@nterm
    nchar <- length(data)
    if (is.null(names(data))) names(data) <- paste("char", 1:nchar, sep = "")
    dat.tmp <- lapply(data, function(x) x[tree@nodes])
    dat <- lapply(dat.tmp, function(y) y[tree@term])
	dat <- dat[[1]]
	nterm <- Ntip(ape.tre)
	nchar <- 1
	w <- matrix(data=1,nrow=nterm,ncol=1)
	n <- length(dat)
	otree.df=as(tree,"data.frame")
	ordering=otree.df[n:nrow(otree.df),'labels']
	b <- vcv.phylo(ape.tre)
	b.match <- match(row.names(b),ordering)
	b <- b[b.match,b.match]
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

split.or.merge <- function(cur.delta, cur.values, phy, theta=FALSE) { 
# mod 09.14.2010 JM Eastman
	bb=cur.delta
	vv=cur.values
	names(vv)<-names(bb)<-phy$edge[,2]
	new.bb=choose.one(bb)
	new.vv=vv
	
	s=names(new.bb[bb!=new.bb])
	all.shifts=as.numeric(names(bb[bb>0]))
	all.D=get.descendants.of.node(s, phy)
	untouchable=unlist(lapply(all.shifts[all.shifts>s], function(x)get.descendants.of.node(x,phy)))
	remain.unchanged=union(all.shifts, untouchable)
	
	marker=match(s, names(new.vv))
	nn=length(vv)
	K=length(table(vv))
	N=Ntip(phy)
	
	if(sum(new.bb)>sum(bb)) {			# add transition: SPLIT
		if(theta) {							# assign new theta to current node
			new.vv[marker] <- nr <- adjust.value(new.vv[marker])
		} else {							# assign new rate to current node
			new.vv[marker] <- nr <- adjust.rate(new.vv[marker])
		}
		if(length(remain.unchanged)==0) {	# assign new value to all descendants
			new.vv[match(all.D,names(new.vv))] = nr
		} else {							# assign new value to all 'open' descendants 
			new.vv[match(all.D[!(all.D%in%remain.unchanged)],names(new.vv))] = nr
		}
		lnHastingsRatio = log((K+1)/(2*N-2-K)) ### from Drummond and Suchard 2010: where N is tips, K is number of local parms in tree
		lnPriorRatio = log(ppois(K+1,log(2),lower.tail=FALSE)/ppois(K,log(2),lower.tail=FALSE))
	} else {							# drop transition: MERGE
		anc = get.ancestor.of.node(s, phy)
		if(!is.root(anc, phy)) {			# substitute ancestral value for current
			new.vv[match(s, names(new.vv))] <- nr <- new.vv[match(anc, names(new.vv))] 
		} else {							
			if(theta) {						# assign new theta to first descendant of root
				new.vv[marker] <- nr <- adjust.value(new.vv[marker])
			} else {						# assign new rate to first descendant of root
				new.vv[marker] <- nr <- adjust.rate(new.vv[marker])
			}		
		}
		if(length(remain.unchanged)==0) {	# assign new value to all descendants
			new.vv[match(all.D, names(new.vv))] = nr
		} else {							# assign new value to all 'open' descendants
			new.vv[match(all.D[!(all.D%in%remain.unchanged)],names(new.vv))] = nr
		}
		lnHastingsRatio = log((2*N-2-K+1)/K) ### from Drummond and Suchard 2010: where N is tips, K is number of local parms in tree
		lnPriorRatio = log(ppois(K-1,log(2),lower.tail=FALSE)/ppois(K,log(2),lower.tail=FALSE))
		
	}
	return(list(new.values=new.vv, new.delta=new.bb, lnHastingsRatio=lnHastingsRatio, lnPriorRatio=lnPriorRatio))
}

assign.root.theta <- function(cur.theta, phy) {
	root=phy$edge[1,1]
	desc=phy$edge[which(phy$edge[,1]==root),2]
	root.theta=cur.theta[sample(which(phy$edge[,2]==desc),1)]
	root.theta
}

get.descendants.of.node <- function(node, phy) {
	descendants.of.node <- function(node, phy) {
		current.node=node[node%in%phy$edge[,1]]
		if (length(current.node)>0) {
			desc=c(sapply(current.node, function(x) {return(phy$edge[which(phy$edge[,1]==x),2])}))
			write(desc,file="TEMP",append=TRUE)
			descendants.of.node(desc,phy)
		}	
	}
	descendants.of.node(node, phy)
	if(file.exists("TEMP")) {
		out=unlist(strsplit(readLines("TEMP")," "))
		unlink("TEMP")
		return(sort(unique(as.numeric(out))))
	}
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

apply.BMM.to.apetree<-function(apetree, rates) {
	apetree$edge.length=apetree$edge.length*rates
	return(apetree)
}

adjust.value <- function(value) {
	vv=value
	rch <- runif(1, min = -0.5, max = 0.5)
	vv <- vv + rch
	vv
}

adjust.rate <- function(rate) {
	rr=rate
	r=log(rr)
	rch=runif(1, min=-1, max=1)
	nr=r+rch
	rr=exp(nr)
	rr
}
	
prepare.ouchdata <- function(phy, data) {
	td <- treedata(phy, data, sort = T)	
	ouch.tre <- ape2ouch(td$phy->ape.tre, scale=FALSE)		
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
	return(list(ouch.tre=ouch.tre, ouch.dat=ouch.dat, ape.tre=td$phy, orig.dat=td$data[,1]))
}

summarize.run<-function(cOK, cProp, cur.model) {
	names(cProp)<-names(cOK)<-c("rates", "rcats", "root", "alpha", "sRegimes", "theta")
	names(cProp)=paste("prop.",names(cProp),sep="")
	names(cOK)=paste("okay.",names(cOK),sep="")
	cat("\n")
	if(cur.model=="BM") {
		print(cProp[1:3])
		print(cOK[1:3])
	} else {
		print(cProp[c(1,2,4,5,6)])
		print(cOK[c(1,2,4,5,6)])
	}
}

generate.log<-function(bundled.parms, cur.model, file, init=FALSE) {
#	bundled.parms=list(gen=i, mrate=mean(cur.rates), cats=max(cur.delta.rates), root=cur.root, alpha=cur.alpha, reg=max(cur.regimes), theta=mean(cur.theta), lnL=curr.lnL)
	res=bundled.parms

	if(init) {
		if(cur.model=="BM") msg=paste("gen", "model", "mean.rate", "rates", "root", "lnL", sep="\t") else msg=paste("gen", "model", "mean.rate", "rates", "alpha", "regimes", "mean.theta", "lnL", sep="\t")
	} else {
		if(cur.model=="BM") {
			msg<-paste(res$gen, sQuote(cur.model), sprintf("%.3f", res$mrate), res$cats, sprintf("%.3f", res$root), sprintf("%.3f", res$lnL),sep="\t")
		} else {
			msg<-paste(res$gen, sQuote(cur.model), sprintf("%.3f", res$mrate), res$cats, sprintf("%.4f", res$alpha), res$reg, sprintf("%.3f", res$theta), sprintf("%.3f", res$lnL),sep="\t")
		}
	}
	write(msg, file=file, append=TRUE)
}

generate.error.message<-function(i, mod.cur, mod.new, lnR, errorLog, init=FALSE) {
	if(init) {
		initmsg = write(paste("gen", "curr.lnL", "new.lnL", "lnR", sep="\t"), file=errorLog)
	} else {
		write(paste(i, sprintf("%.3f", mod.cur$lnL), sprintf("%.3f", mod.new$lnL), sprintf("%.3f", lnR), sep="\t"), file=errorLog)
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
	