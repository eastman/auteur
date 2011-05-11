#rjmcmc proposal mechanism: proposal mechanism for traversing model complexity (by one parameter at a time)
#author: JM EASTMAN 2010

splitormerge <-
function(cur.delta, cur.values, phy, node.des, lambda=lambda, logspace=TRUE, internal.only=FALSE) { 
	bb=cur.delta
	vv=cur.values
	names(vv)<-names(bb)<-phy$edge[,2]
	new.bb=choose.one(bb, phy=phy, internal.only)
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
	K=sum(bb)
	N=Ntip(phy)
	
	ncat=sum(bb)
	cur.vv=as.numeric(vv[marker])
	ca.vv=length(which(vv==cur.vv))
	
	if(sum(new.bb)>sum(bb)) {			# add transition: SPLIT
		decision="split"
		n.desc=sum(!all.D%in%remain.unchanged)+1
		n.split=sum(vv==cur.vv)-n.desc
		if(!logspace) {
			u=splitvalue(cur.vv=cur.vv, n.desc=n.desc, n.split=n.split, factor=lambda) 
		} else {
			u=splitrate(cur.vv=cur.vv, n.desc=n.desc, n.split=n.split, factor=lambda)
		}
		nr.split=u$nr.split
		nr.desc=u$nr.desc
		new.vv[vv==cur.vv]=nr.split
		if(length(remain.unchanged)==0) {	# assign new value to all descendants
			new.vv[match(all.D,names(new.vv))] = nr.desc
		} else {							# assign new value to all 'open' descendants 
			new.vv[match(all.D[!(all.D%in%remain.unchanged)],names(new.vv))] = nr.desc
		}
		new.vv[match(s, names(new.vv))]=nr.desc
		lnHastingsRatio = log((K+1)/(2*N-2-K)) ### from Drummond and Suchard 2010: where N is tips, K is number of local parms in tree
#		lnPriorRatio = log(ptpois(K+1,lambda,nn)/ptpois(K,lambda,nn))
		lnPriorRatio = dpois(K+1,lambda,log=TRUE)-dpois(K,lambda,log=TRUE)
		
	} else {							# drop transition: MERGE
		decision="merge"
		anc = get.ancestor.of.node(s, phy)
		if(!is.root(anc, phy)) {			# base new rate on ancestral rate of selected branch
			anc.vv=as.numeric(vv[match(anc,names(vv))])
			na.vv=length(which(vv==anc.vv))
			nr=(anc.vv*na.vv+cur.vv*ca.vv)/(ca.vv+na.vv)
			new.vv[vv==cur.vv | vv==anc.vv]=nr
		} else {							# if ancestor of selected node is root, base new rate on sister node
			sister.tmp=get.desc.of.node(anc,phy)
			sister=sister.tmp[sister.tmp!=s]
			sis.vv=as.numeric(vv[match(sister,names(vv))])
			ns.vv=length(which(vv==sis.vv))
			nr=(sis.vv*ns.vv+cur.vv*ca.vv)/(ca.vv+ns.vv)
			new.vv[vv==cur.vv | vv==sis.vv]=nr			
		}
		lnHastingsRatio = log((2*N-2-K+1)/K) ### from Drummond and Suchard 2010: where N is tips, K is number of local parms in tree
#		lnPriorRatio = log(ptpois(K-1,lambda,nn)/ptpois(K,lambda,nn))
		lnPriorRatio = dpois(K-1,lambda,log=TRUE)-dpois(K,lambda,log=TRUE)

		
	}
	
	return(list(new.delta=new.bb, new.values=new.vv, lnHastingsRatio=lnHastingsRatio, lnPriorRatio=lnPriorRatio, decision=decision))
}


#general phylogenetic utility for selecting one element of a list of values, each of which can be associated with the edges of a phylogeny (phy$edge[,2])
#author: JM EASTMAN 2011

choose.one <-
function(cur.delta, phy, internal.only=FALSE)
# updated 02.26.2011 to weight by branch length (rather than equal weight for each edge)
{
	edge.prob=0.75
	e=phy$edge.length
	ee=cumsum(e/sum(e))
	bb=cur.delta
	if(internal.only) {
		names(bb)=phy$edge[,2]
		while(1) {
			if(runif(1)<edge.prob) {
				s=min(which(runif(1)<ee))
			} else {
				s=sample(1:length(bb),1)
			}
			if(as.numeric(names(bb[s]))>Ntip(phy)) break()
		}
	} else {
		if(runif(1)<edge.prob) {
			s=min(which(runif(1)<ee))
		} else {
			s=sample(1:length(bb),1)
		}
	}
	bb[s]=1-bb[s]
	bb
}


#rjmcmc proposal mechanism: split a rate into two while maintaining the mean
#author: JM EASTMAN 2010

splitrate <-
function(cur.vv, n.desc, n.split, factor=log(2)){
	dev=cur.vv-exp(log(cur.vv)+runif(1, -factor, factor))
	nr.desc=cur.vv+dev/n.desc
	nr.split=cur.vv-dev/n.split
	return(list(nr.desc=nr.desc, nr.split=nr.split))
}


#rjmcmc proposal mechanism: split a value into two while maintaining the mean
#author: JM EASTMAN 2010

splitvalue <-
function(cur.vv, n.desc, n.split, factor=log(2)){
	dev=cur.vv-adjustvalue(cur.vv, factor)
	nr.desc=cur.vv + dev/n.desc
	nr.split=cur.vv - dev/n.split
	return(list(nr.desc=nr.desc, nr.split=nr.split))
}


#rjmcmc proposal mechanism for updating lineage-specific relative rates (i.e., a subvector swap)
#author: JM EASTMAN 2010

tune.rate <-
function(rates, prop.width) {
	ss=sample(rates, 1)
	ww=which(rates==ss)
	nn=adjustrate(ss,prop.width)
	rates[ww]=nn
	return(rates)
}


#rjmcmc proposal mechanism for updating a single numeric class from a vector
#author: JM EASTMAN 2010

tune.value <-
function(values, prop.width) {
	ss=sample(values, 1)
	ww=which(values==ss)
	nn=adjustvalue(ss,prop.width)
#	nn=proposal.slidingwindow(value, prop.width, min=NULL, max=NULL)
	values[ww]=nn
	return(values)
}


#scales VCV matrix of phylogeny by relative evolutionary rates under Brownian motion
#author: JM EASTMAN 2010

updatevcv <-
function(ape.tre, new.rates) {
	n=ape.tre
	n$edge.length=n$edge.length*new.rates
	vv=vmat(n)
	return(vv)
}



#rjmcmc proposal mechanism: scale a value given a proposal width ('prop.width')
#author: JM EASTMAN 2010

adjustvalue <-
function(value, prop.width) {
# mod 10.20.2010 JM Eastman
## FIXME: defunct

	vv=value
	if(runif(1)<0.95 | length(unique(vv))==1) {
		rch <- runif(1, min = -prop.width, max = prop.width)
		return(vv[]+rch)
	} else {
		return((vv-mean(vv))*runif(1,min=-prop.width, max=prop.width)+mean(vv))
	}
}


#rjmcmc proposal mechanism: scale a rate (limited to only positive values) given a proposal width ('prop.width')
#author: JM EASTMAN 2010
## FIXME: defunct

adjustrate <-
function(rate, prop.width) {
# mod 05.06.2011 JM Eastman
	vv=rate
	v=log(vv)
	rch <- runif(1, min = -prop.width, max = prop.width)
	v=exp(v[]+rch)
	v
}


#mcmc prior calculation: compute ln-prior ratio for exponential-distributed values
#author: JM EASTMAN 2011


priorratio.exp <- function(cur.values, cur.delta, new.values, new.delta, mean){
	
	proc.priorratio.exp=function(v,d,m){
		t=v[d==1]
		p=ifelse(length(t), sum(dexp(t, rate=1/m, log=TRUE)), 1)
		p
	}
	
	lnPriorRatio = proc.priorratio.exp(new.values, new.delta, mean) - proc.priorratio.exp(cur.values, cur.delta, mean)
	return(lnPriorRatio)
}

#rjmcmc proposal mechanism: adjust a value (within bounds)
#author: JM EASTMAN 2011

proposal.slidingwindow <- function(value, prop.width, min=NULL, max=NULL){
	u=runif(1)
	v=value+(u-0.5)*prop.width
	if((!is.null(min) | !is.null(max))){
		
		# reflect if out-of-bounds
		if(any(v>max)) {
			v[v>max]=max-(v[v>max]-max)
		} else if(any(v<min)){
			v[v<min]=min-(v[v<min]-min)
		}
		
	}	
	return(list(v=v, lnHastingsRatio=0))
}


#rjmcmc proposal mechanism: scale a value with asymmetrically drawn multiplier
#author: JM EASTMAN 2011

proposal.multiplier <- function(value, prop.width){
	tmp=c(prop.width, 1/prop.width)
	a=min(tmp)
	b=max(tmp)
	lambda=2*log(b)
	u=runif(1)
	m=exp(lambda*(u-0.5))
	return(list(v=value*m, lnHastingsRatio=log(u)))
}
