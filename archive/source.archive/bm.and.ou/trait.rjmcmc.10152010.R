## Written by LJ HARMON & AH: 2008
## Modified by JM EASTMAN: 09.2010

rjMCMC.trait<-function (phy, data, ngen=1000, sampleFreq=100, model=c("OU","BM"), probAlpha=0.05, probMergeSplit=0.05, probRoot=0.1, theta.rates.ratio=1, heat=1, fileBase="result", simplestart=FALSE) 
{
	if(length(model)>1)stop(paste("Please specify either ", sQuote("OU"), " or ", sQuote("BM")))
## MORE ERROR CHECKING NEEDED HERE ## 
	require(ouch)
	require(geiger)
	
### prepare objects for rjMCMC
	cur.model		<- model								

	dataList		<- prepare.ouchdata(phy, data)			
	ape.tre			<- cur.tre <- dataList$ape.tre			
	orig.dat		<- dataList$orig.dat						
	ouch.tre		<- dataList$ouch.tre					
	ouch.dat		<- dataList$ouch.dat						
	nn				<- length(ape.tre$edge.length)	
	node.des		<- sapply(unique(c(ape.tre$edge[1,1],ape.tre$edge[,2])), function(x) get.descendants.of.node(x, ape.tre))
	names(node.des) <- c(ape.tre$edge[1,1], unique(ape.tre$edge[,2]))
	ll <- list()										

# initialize parameters
	if(simplestart) {
		init.rate	<- list(values=rep(fitContinuous(phy,data)$Trait1$beta,length(phy$edge.length)),delta=rep(0,length(phy$edge.length)))
		init.theta	<- list(values=rep(mean(data),length(phy$edge.length)), delta=rep(0,length(phy$edge.length)))
	} else {
		init.rate	<- generate.starting.point(data,ape.tre,node.des,theta=FALSE)
		init.theta	<- generate.starting.point(data,ape.tre,node.des,theta=TRUE)

	}
	
    cur.rates		<- init.rate$values			
	cur.delta.rates	<- init.rate$delta
	cur.root		<- adjust.value(mean(orig.dat))
	cur.alpha		<- adjust.rate(0.001)
	cur.theta		<- init.theta$values							
	cur.delta.theta	<- init.theta$delta	
	cur.root.theta	<- assign.root.theta(cur.theta,ape.tre)
	cur.regimes		<- as.factor(c(cur.root.theta,cur.theta)) 
	cur.tre			<- apply.BMM.to.apetree(ape.tre, cur.rates)
	
	if(model=="OU")	{
		mod.cur = new.ou.lik.fn(cur.rates, cur.alpha, cur.regimes, cur.theta, cur.tre, ouch.dat)
	} else {
		mod.cur = new.bm.lik.fn(cur.rates, cur.root, cur.tre, orig.dat)
	}


# proposal counts
	nRateProp =		0	
	nRateSwapProp = 0									
	nRcatProp =		0										
	nRootProp =		0
	nAlphaProp =	0										
	nTcatProp =		0										
	nThetaProp =	0
	nThetaSwapProp= 0
	nRateOK =		0											
	nRateSwapOK =	0
	nRcatOK =		0											
	nRootOK =		0
	nAlphaOK =		0										
	nTcatOK =		0
	nThetaOK =		0
	nThetaSwapOK =	0

	cOK=c(												
		nRateOK,
		nRateSwapOK,
		nRcatOK, 
		nRootOK, 
		nAlphaOK,
		nThetaOK,
		nThetaSwapOK,
		nTcatOK 
	)		
	
	tickerFreq=ceiling(ngen/30)
	
# file handling
	errorLog=file(paste(getwd(),paste(cur.model, fileBase, "rjTraitErrors.log",sep="."),sep="/"),open='w+')
	generate.error.message(i=NULL, mod.cur=NULL, mod.new=NULL, lnR=NULL, errorLog=errorLog, init=TRUE)
	runLog=file(paste(getwd(),paste(cur.model, fileBase, "rjTrait.log",sep="."),sep="/"),open='w+')
	generate.log(bundled.parms=NULL, cur.model, file=runLog, init=TRUE)

### Begin rjMCMC
    for (i in 1:ngen) {
        lnLikelihoodRatio <- lnHastingsRatio <- lnPriorRatio <- 0
		startparms = c(nRateProp, nRateSwapProp, nRcatProp, nRootProp, nAlphaProp, nThetaProp, nThetaSwapProp, nTcatProp)
	
## OU IMPLEMENTATION ## 
		if(cur.model=="OU") {
			if (runif(1) < (2 * probMergeSplit)) {
				if(runif(1)<(theta.rates.ratio/(theta.rates.ratio+1))) {	# adjust theta categories
					nr=split.or.merge(cur.delta.theta, cur.theta, ape.tre, node.des, theta=TRUE)
					new.alpha=cur.alpha
					new.theta=nr$new.values ##
					new.delta.theta=nr$new.delta ##
					new.root.theta<-nrt<-assign.root.theta(new.theta, ape.tre) ## 
					new.regimes=as.factor(c(nrt,new.theta)) ##
					new.rates=cur.rates 
					new.delta.rates=cur.delta.rates 
					nTcatProp=nTcatProp+1
				} else {							# adjust rate categories
					nr=split.or.merge(cur.delta.rates, cur.rates, ape.tre, node.des, theta=FALSE)
					new.alpha=cur.alpha
					new.root.theta=cur.root.theta
					new.theta=cur.theta
					new.delta.theta=cur.delta.theta
					new.regimes=cur.regimes
					new.rates=nr$new.values ##
					new.delta.rates=nr$new.delta ##
					nRcatProp=nRcatProp+1
				}
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
			} else { # neither split nor merge
				if(runif(1)<probAlpha) {			# adjust alpha
					new.alpha=adjust.rate(cur.alpha) ##
					new.root.theta=cur.root.theta
					new.theta=cur.theta
					new.delta.theta=cur.delta.theta
					new.regimes=cur.regimes 
					new.rates=cur.rates
					new.delta.rates=cur.delta.rates
					nAlphaProp=nAlphaProp+1
				} else { 
					if(runif(1)>(theta.rates.ratio/(theta.rates.ratio+1))) {
						if(runif(1)<0.2) {
							new.alpha=cur.alpha		# adjust rates
							new.root.theta=cur.root.theta
							new.theta=cur.theta
							new.delta.theta=cur.delta.theta
							new.regimes=cur.regimes 
							new.rates=adjust.rate(cur.rates) ##
							new.delta.rates=cur.delta.rates
							nRateProp=nRateProp+1
						} else if(sum(cur.delta.rates)!=0) {					# swap rates
							nr=switcheroo(cur.delta.rates, cur.rates, ape.tre, node.des, theta=FALSE)
							new.alpha=cur.alpha		
							new.root.theta=cur.root.theta
							new.theta=cur.theta
							new.delta.theta=cur.delta.theta
							new.regimes=cur.regimes 
							new.rates=nr$new.values ##
							new.delta.rates=nr$new.delta ##
							nRateSwapProp=nRateSwapProp+1
						}
					} else {	
						if(runif(1)<0.2 | sum(cur.delta.theta)==0) {			# adjust theta
							new.alpha=cur.alpha 
							new.theta=adjust.value(cur.theta) ##
							new.delta.theta=cur.delta.theta
							new.root.theta<-nrt<-assign.root.theta(new.theta, ape.tre) ## 
							new.regimes=as.factor(c(nrt,new.theta)) ##
							new.rates=cur.rates
							new.delta.rates=cur.delta.rates
							nThetaProp=nThetaProp+1	
						} else {					# swap theta
							nr=switcheroo(cur.delta.theta, cur.theta, ape.tre, node.des, theta=TRUE)
							new.alpha=cur.alpha 
							new.theta=nr$new.values ##
							new.delta.theta=nr$new.delta ##
							new.root.theta<-nrt<-assign.root.theta(new.theta, ape.tre) ## 
							new.regimes=as.factor(c(nrt,new.theta)) ##
							new.rates=cur.rates
							new.delta.rates=cur.delta.rates
							nThetaSwapProp=nThetaSwapProp+1	
						}
					}
				}
			}
			
			new.tre=apply.BMM.to.apetree(ape.tre, new.rates)
			mod.new=new.ou.lik.fn(new.rates, new.alpha, new.regimes, new.theta, new.tre, ouch.dat)
		} else if(cur.model=="BM"){
## BM IMPLEMENTATION ##
			if (runif(1) < (2 * probMergeSplit)) {							# adjust rate categories
				nr=split.or.merge(cur.delta.rates, cur.rates, ape.tre, node.des, theta=FALSE)
				new.rates=nr$new.values
				new.delta.rates=nr$new.delta
				new.root=cur.root
				nRcatProp=nRcatProp+1
				lnHastingsRatio=nr$lnHastingsRatio
				lnPriorRatio=nr$lnPriorRatio
			} else { 
				if(runif(1)<probRoot) {										# adjust root
					new.root=adjust.value(cur.root)
					new.rates=cur.rates
					new.delta.rates=cur.delta.rates
					nRootProp=nRootProp+1
				} else {													
					if(runif(1)>0.1 & sum(cur.delta.rates)!=0) {			# swap rate classes
						nr=switcheroo(cur.delta.rates, cur.rates, ape.tre, node.des, theta=FALSE)
						new.rates=nr$new.values ##
						new.delta.rates=nr$new.delta ##
						new.root=cur.root
						nRateSwapProp=nRateSwapProp+1
					} else {												# adjust rates
						new.rates=adjust.rate(cur.rates)
						new.delta.rates=cur.delta.rates
						new.root=cur.root
						nRateProp=nRateProp+1
					} 
				}
			}
			
			new.tre=apply.BMM.to.apetree(ape.tre, new.rates)
			mod.new=try(new.bm.lik.fn(new.rates, new.root, new.tre, orig.dat),silent=TRUE)
		}
		if(inherits(mod.new, "try-error")) {mod.new=as.list(mod.new); mod.new$lnL=Inf}
		if(!is.infinite(mod.new$lnL)) {
			lnLikelihoodRatio = mod.new$lnL - mod.cur$lnL
		} else {
			lnLikelihoodRatio = -Inf
			failures=paste("./failed", model, "parameters/",sep=".")
			failed.trees=paste(failures, "trees", sep="/")
			lapply(list(failures, failed.trees), function(x) if(!file.exists(x)) dir.create(x))
				
			pdf(file=paste(failed.trees, paste(i,"failed.model.pdf",sep="."),sep="/"))
				plot(new.tre)
			dev.off()
			
			write.table(matrix(c(i, new.rates),nrow=1),file=paste(failures,"failed.rates.txt",sep="/"),append=TRUE,col=FALSE,row=FALSE,quote=FALSE)
			write.table(matrix(c(i, new.root),nrow=1),file=paste(failures,"failed.root.txt",sep="/"),append=TRUE,col=FALSE,row=FALSE,quote=FALSE)
		}
		
	# compare likelihoods
		endparms = c(nRateProp, nRateSwapProp, nRcatProp, nRootProp, nAlphaProp, nThetaProp, nThetaSwapProp, nTcatProp)
		
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
			mod.cur <- mod.new
			cOK <- determine.accepted.proposal(startparms, endparms, cOK)
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
			bundled.parms=list(gen=i, mrate=exp(mean(log(cur.rates))), cats=sum(cur.delta.rates)+1, root=cur.root, alpha=cur.alpha, reg=(sum(cur.delta.theta)+1), theta=mean(cur.theta), lnL=curr.lnL)
			generate.log(bundled.parms, cur.model, file=runLog)
			ll$lnL=rbind(ll$lnL, curr.lnL)
			ll$rates=rbind(ll$rates, cur.rates) 
			if(cur.model=="OU") {
				ll$theta=rbind(ll$theta, cur.theta)
				ll$alpha=rbind(ll$alpha, cur.alpha)
				ll$root=rbind(ll$root, cur.root.theta)
			}
			if(cur.model=="BM") {
				ll$root=rbind(ll$root, cur.root)
			}
			if(i==ngen) {
				if(cur.model=="OU") dimnames(ll$theta)[[2]]=ape.tre$edge[,2]
				dimnames(ll$rates)[[2]]=ape.tre$edge[,2]
			}
		}
	}
	
# End rjMCMC
	close(errorLog)
	close(runLog)
	summarize.run(cOK, endparms, cur.model)	
	return(list(results=ll,data=orig.dat,tree=ape.tre))
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

new.bm.lik.fn <- function(rates, root, ape.tre, orig.dat) { # from OUCH 2.7-1:::brown()	
# mod 10.12.2010 JM Eastman
# returning INF det(b) on some of these trees -- especially large trees
	tree=ape.tre
	y=orig.dat
	
	b <- vcv.phylo(ape.tre)
	w <- rep(root, nrow(b))
	num <- -t(y-w)%*%solve(b)%*%(y-w)/2
	den <- sqrt((2*pi)^length(y)*det(b))
	list(
		 root=root,
		 lnL = (num-log(den))[1,1]
		 )	
}

split.or.merge <- function(cur.delta, cur.values, phy, node.des, theta=FALSE) { 
# mod 09.17.2010 JM Eastman
	bb=cur.delta
	vv=cur.values
	names(vv)<-names(bb)<-phy$edge[,2]
	new.bb=choose.one(bb)
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
	
	if(sum(new.bb)>sum(bb)) {			# add transition: SPLIT
#		print("s")
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
#		print("m")
		anc = get.ancestor.of.node(s, phy)
		if(!is.root(anc, phy)) {			# substitute ancestral value for current
			new.vv[match(s, names(new.vv))] <- nr <- new.vv[match(anc, names(new.vv))] 
		} else {
			sister.tmp=get.desc.of.node(anc,phy)
			sister=sister.tmp[sister.tmp!=s]
			if(bb[match(sister, names(bb))]==0) {
				nr=as.numeric(vv[match(sister, names(vv))])
			} else {
				if(theta) {						# assign new theta to first descendant of root
					new.vv[marker] <- nr <- adjust.value(new.vv[marker])
				} else {						# assign new rate to first descendant of root
					new.vv[marker] <- nr <- adjust.rate(new.vv[marker])
				}
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
	if(((length(table(new.vv))-1)-sum(new.bb))!=0) fail=TRUE else fail=FALSE
	return(list(new.delta=new.bb, new.values=new.vv, lnHastingsRatio=lnHastingsRatio, lnPriorRatio=lnPriorRatio, fail=fail))
}

switcheroo <- function(cur.delta, cur.values, phy, node.des, theta=FALSE) { 
	bb=cur.delta
	vv=cur.values
	new.bb=sample(bb)
	names(vv)<-names(bb)<-names(new.bb)<-phy$edge[,2]
	new.vv=vv
	
	old.shifts=as.numeric(names(bb[bb!=0]))
	old.values=vv[match(old.shifts, names(vv))]
	
	anc=lapply(old.shifts, function(x)get.ancestor.of.node(x,phy))
	if(!is.root(min(unlist(anc)),phy)) {
	# find most ancestral value in tree
		anc.value=unname(vv[match(min(unlist(anc)),names(vv))]) 
	} else {
		if(theta) anc.value=mean(sapply(old.values, function(x) adjust.value(x))) else anc.value=mean(sapply(old.values, function(x) adjust.rate(x)))
	}
	all.desc.anc=node.des[[which(names(node.des)==min(unlist(anc)))]]
	new.vv[match(all.desc.anc,names(new.vv))] = anc.value

	new.shifts=sort(as.numeric(names(new.bb[new.bb!=0])))
	new.tips=new.shifts[new.shifts<=Ntip(phy)]
	new.ints=new.shifts[new.shifts>Ntip(phy)]
	new.shifts=c(new.ints, new.tips)
	
	new.values=as.numeric(old.values)[sample(1:length(old.values))]
	for(z in 1:length(new.shifts)) {
		node=new.shifts[z]
		new.vv[match(node, names(new.vv))]=new.values[z]
		new.vv[match(node.des[[which(names(node.des)==node)]],names(new.vv))]=new.values[z]
	} 
	return(list(new.delta=new.bb, new.values=new.vv))
}

generate.starting.point <- function(data, phy, node.des, theta=FALSE) { 
	nn=length(phy$edge.length)
	ntip=Ntip(phy)
	nshifts=rtpois(1,log(ntip),nn)
	if(nshifts!=0) bb=sample(c(rep(1,nshifts),rep(0,nn-nshifts)),replace=FALSE) else bb=rep(0,nn)
	names(bb)=phy$edge[,2]
	shifts=as.numeric(names(bb)[bb==1])
	if(theta) {
		min.max=c(adjust.value(min(data)), adjust.value(max(data)))
		values=runif(sum(bb)+1, min=min.max[min.max==min(min.max)], max=min.max[min.max==max(min.max)])
	} else {
		init.rate=fitContinuous(phy,data)[[1]]$beta
		min.max=c(adjust.rate(init.rate), adjust.rate(init.rate))
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
#		d=get.descendants.of.node(internal.shifts[i],phy)
				vv[match(c(internal.shifts[i], d), names(vv))]=values[i]
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

get.descendants.of.node <- function(node, phy) {
	storage=list()
	if(node%in%phy$edge[,1]) {
		desc=c(sapply(node, function(x) {return(phy$edge[which(phy$edge[,1]==x),2])}))
		while(any(desc%in%phy$edge[,1])) {
			desc.true=desc[which(desc%in%phy$edge[,1])]
			desc.desc=lapply(desc.true, function(x) return(get.desc.of.node(x, phy)))
			
			storage=list(storage,union(desc,desc.desc))
			
			desc=unlist(desc.desc)
		}
		storage=sort(unique(unlist(union(desc,storage))))
	}
	return(storage)
}

get.desc.of.node <- function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
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

adjust.value <- function(value) {
	vv=value
	vv.levels=levels(factor(vv))
	rch <- runif(length(vv.levels), min = -1, max = 1)
	for(i in 1:length(vv.levels)){
		vv[which(vv==vv.levels[i])]=vv[which(vv==vv.levels[i])]+rch[i]
	}
	vv
}

adjust.rate <- function(rate) {
	vv=rate
	v=log(vv)
	vv.levels=levels(factor(vv))
	rch <- runif(length(vv.levels), min = -1, max = 1)
	for(i in 1:length(vv.levels)){
		v[which(vv==vv.levels[i])]=v[which(vv==vv.levels[i])]+rch[i]
	}
	v=exp(v)
	v
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

summarize.run<-function(cOK, endparms, cur.model) {
#	end = c(nRateProp, nRateSwapProp, nRcatProp, nRootProp, nAlphaProp, nThetaProp, nThetaSwapProp, nTcatProp)

	names(endparms)<-names(cOK)<-c("rates", "rate.swap", "rate.catg", "root", "alpha", "theta", "theta.swap", "theta.catg")
	names(endparms)=paste("pr.",names(endparms),sep="")
	names(cOK)=paste("ok.",names(cOK),sep="")
	cat("\n")
	if(cur.model=="BM") {
		print(endparms[paste("pr.", c("rates", "rate.swap", "rate.catg", "root"),sep="")])
		print(cOK[paste("ok.", c("rates", "rate.swap", "rate.catg", "root"),sep="")])
	} else {
		print(endparms[paste("pr.", c("rates", "rate.swap", "rate.catg", "alpha", "theta", "theta.swap", "theta.catg"),sep="")])
		print(cOK[paste("ok.", c("rates", "rate.swap", "rate.catg", "alpha", "theta", "theta.swap", "theta.catg"),sep="")])
	}
}

generate.log<-function(bundled.parms, cur.model, file, init=FALSE) {
#	bundled.parms=list(gen=i, mrate=mean(cur.rates), cats=sum(cur.delta.rates)+1, root=cur.root, alpha=cur.alpha, reg=sum(cur.delta.theta)+1, theta=mean(cur.theta), lnL=curr.lnL)
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
				 

				 
	