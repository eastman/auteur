
# MAIN FUNCTION for LEVY PROCESS
levyevolver=function(t, alpha, sigmasq.brown, sigma.jump, lambda.jump){
	y<-B<-rnorm(1, mean=alpha, sd=sqrt(sigmasq.brown*t))
	jumps=rpois(1, lambda.jump*t)
	if(jumps) {
		y=y+(L<-sum(rnorm(jumps, mean=y, sd=sigma.jump)))
	} else {
		L=alpha
	}
	return(list(y=y, Brownian.effect=B-alpha, Levy.effect=L-alpha))
}

phenogram<-function(phy,model=c("levy","bm"),alpha=0,rate=0.01,lambda=0.5,sigma=0.5,increments=100,plot=TRUE,node.cex=2,...){
# written by JM Eastman 02.22.2011 as an extension of fastBM() by LJ Revell, under a Levy process (as a compound Brownian-Poisson process)
	
# plotting
	add.transparency <-
	function (col, alpha) {
		tmp <- col2rgb(col)/255
		rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
	}
	
# find segments of time within a vector 
	get.lengths=function(yy) {
		yy[2:length(yy)]-yy[1:(length(yy)-1)]
	}

	
	tree=phy 
	bb=branching.times(tree)
	abstimes=max(bb)-bb
	if(!is.null(increments)) increment=max(bb)/increments else increment=NULL
	
	n<-length(tree$tip.label)
	
## first simulate changes along each branch
	if(model=="levy") {
		## ENSURE THAT DRAW FROM LEVY PROCESS IS ACCURATELY MODELLED ##
				# evolves trait by Brownian motion as usual
				# adds additional deviate according to Poisson-distributed jumps and a given jump size
		x<-lapply(1:length(tree$edge.length), function(x) {
				  if(is.numeric(increment)) y=seq(0,tree$edge.length[x],by=increment) else y=0
				  if(length(y)==1) {
#					return(sample(c(-1,1),1)*rlevy(n=1, m=0, s=sqrt(rate*tree$edge.length[x])))
					return(levyevolver(t=tree$edge.length[x], alpha=alpha, sigmasq.brown=rate, sigma.jump=sigma, lambda.jump=lambda)$y)
				  } else {
#					return(sapply(get.lengths(y), function(z) return(sample(c(-1,1),1)*rlevy(n=1, m=0, s=sqrt(rate*z)))))
					return(sapply(get.lengths(y), function(z) return(levyevolver(t=z, alpha=alpha, sigmasq.brown=rate, sigma.jump=sigma, lambda.jump=lambda)$y)))

				  }
			}
		)
		plot.title="Levy process"
	} else if(model=="bm") {
		x<-lapply(1:length(tree$edge.length), function(x) {
				  if(is.numeric(increment)) y=seq(0,tree$edge.length[x],by=increment) else y=0
				  if(length(y)==1) {
					return(rnorm(1, mean=0, sd=sqrt(rate*tree$edge.length[x])))
				  } else {
					return(sapply(get.lengths(y), function(z) rnorm(n=1,mean=0,sd=sqrt(rate*z))))
				  }
			}
		)		
		plot.title="Brownian"
	}
	
	node.phenotypes=matrix(cbind(tree$edge,0,max(bb)),ncol=4)
	for(i in 1:length(abstimes)) {
		if(names(abstimes)[i]%in%node.phenotypes[,2]) node.phenotypes[which(node.phenotypes[,2]==names(abstimes)[i]),4]=abstimes[i]
	}
	
## accumulate change along each path
	traithistory=list()
	times=list()
	y<-matrix(0,nrow(tree$edge),ncol(tree$edge)+1)
	for(i in 1:length(x)){
		yy=c() # phenotypes
		ntt=c() # node times
		tt=seq(0,tree$edge.length[i],length=length(x[[i]])+1)[-1] # each time increment along a branch
		if(length(tt)==0) tt=tree$edge.length[i] # time where branch is smaller than increment
		if(node.phenotypes[i,1]==n+1) { # ancestor
			start=alpha 
			start.time=0
		} else {						# non-ancestor
			start=node.phenotypes[which(node.phenotypes[,2]==node.phenotypes[i,1]),3]
			start.time=node.phenotypes[which(node.phenotypes[,2]==node.phenotypes[i,1]),4]
		}
		for(j in 1:length(x[[i]])) {
			ntt[j]=start.time+tt[j]
			if(j==1) {
				yy[j]=start+x[[i]][j] 
			} else {
				yy[j]=yy[j-1]+x[[i]][j]
			}
			if(j==length(x[[i]])) node.phenotypes[i,3]=yy[length(yy)]
		}
		traithistory[[i]]=c(start,yy)
		times[[i]]=c(start.time,ntt)
	}
	
	mm=max(abs(alpha-unlist(traithistory)))
	
	# deal with non-ultrametric trees
	for(i in 1:nrow(node.phenotypes)) {
		tip=node.phenotypes[i,2]
		if(!tip%in%tree$edge[,1]) {
			last.edge=node.phenotypes[match(node.phenotypes[i,1],node.phenotypes[,2]),4]
			if(is.na(last.edge)) last.edge=0
			node.phenotypes[i,4]=last.edge+tree$edge.length[which(tree$edge[,2]==tip)]
			
		}
	}
	
	
	
	## PLOTTING ##
	if(plot) {
		plot(x=NULL, y=NULL, xlim=c(0,max(node.phenotypes[,4])), ylim=c(-2*mm+alpha, 2*mm+alpha), bty="n", xlab="time", ylab="phenotypic value",main=paste("evolution by",plot.title,sep=" "))
		for(i in 1:length(x)) {
			lines(times[[i]],traithistory[[i]],col=add.transparency("gray25",0.75),...)		
		}
		for(i in 1:length(x)) {
			points(node.phenotypes[i,4],node.phenotypes[i,3],bg=add.transparency("white",0.75),pch=21,cex=node.cex)
			
		}
		points(0,alpha,bg=add.transparency("white",0.75),pch=21,cex=node.cex)		
	}
	pp=data.frame(node.phenotypes)
	names(pp)=c("ancestor","descendant","phenotype","time")
	
	
	tt=which(pp$descendant<=Ntip(phy))
	tips.dat=pp$phenotype[tt]
	names(tips.dat)=phy$tip.label[pp$descendant[tt]]
	
	return(list(dat=tips.dat, phy=tree, phenotypes=pp))
}
