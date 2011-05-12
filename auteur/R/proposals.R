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
	
	if(sum(new.bb)>sum(bb)) {			## add transition: SPLIT
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
		lnPriorRatio = dpois(K+1,lambda,log=TRUE)-dpois(K,lambda,log=TRUE)
		
	} else {							## drop transition: MERGE
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
		lnPriorRatio = dpois(K-1,lambda,log=TRUE)-dpois(K,lambda,log=TRUE)
	}
	
	return(list(new.delta=new.bb, new.values=new.vv, lnHastingsRatio=lnHastingsRatio, lnPriorRatio=lnPriorRatio, decision=decision))
}



#rjmcmc proposal mechanism: locally select an edge based on subtended or subtending branch lengths
#author: JM EASTMAN 2011

tree.slide <-
function(phy, node, up=NA, use.edges=FALSE){
	if(is.na(up)) up=as.logical(round(runif(1)))
	root=Ntip(phy)+1
	if(up){
		choice.tmp=get.ancestor.of.node(node,phy)
		if(choice.tmp==root){	# choice will need to traverse around root
			sisters.tmp=get.desc.of.node(choice.tmp,phy)
			sisters=sisters.tmp[which(sisters.tmp!=node)]
			if(length(sisters)==1) {
				choice=sisters
				direction="root"
			} else {
				sister.edges=phy$edge.length[match(sisters, phy$edge[,2])]
				if(!use.edges) sister.edges=rep(1,length(sister.edges))
				see=cumsum(sister.edges)
				see=see/max(see)
				choice=sisters[min(which(runif(1)<see))]
				direction="root"
			}
		} else {				# choice is rootward
			choice=choice.tmp
			direction="up"
		}
	} else {					# choice is tipward
		descendants=get.desc.of.node(node,phy)
		if(!length(descendants)) {		# -- choice is a tip
			if(runif(1)<0.5) {
				choice=node
				direction=NULL
			} else {
				choice=tree.slide(phy, node, up=TRUE)$node
				direction="up"
			}
		} else {						# -- choice is an internal node
			desc.edges=phy$edge.length[match(descendants, phy$edge[,2])]
			if(!use.edges) desc.edges=rep(1,length(desc.edges))
			dee=cumsum(desc.edges)
			dee=dee/max(dee)
			choice=descendants[min(which(runif(1)<dee))]
			direction="down"
		}
	}
	
	return(list(node=choice, direction=direction))
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



#proposal utility: move shift-point in tree 
#author: JM EASTMAN 2011

adjustshift <-
function(phy, cur.delta, cur.values){
	delta=cur.delta
	values=cur.values
	
	root=Ntip(phy)+1
	root.des=get.desc.of.node(root, phy)
	names(delta)<-names(values)<-phy$edge[,2]
	shifts=as.numeric(names(which(delta==1)))
	node=shifts[sample(1:length(shifts),1)]
	val=as.numeric(values[which(names(values)==node)])
	
	# select new node
	tmp=tree.slide(phy, node, up=NA)
	newnode=tmp$node
	direction=tmp$direction
	
	# store all current shifts
	shifts=shifts[shifts!=node]
	
#	print(newnode)

	# determine if new shift is non-permissible
	if(newnode%in%c(node,shifts) | all(root.des%in%c(shifts,newnode))) {
#		print("skip")
#		plot(phy, edge.color=ifelse(values==1, 1, ifelse(values==2, "red", "green")))
#		nodelabels(frame="cir",cex=0.75,bg=NA)
		return(list(new.delta=delta, new.values=values, lnHastingsRatio=0))
	} else {		
		
		# update delta
		new.delta=delta
		new.delta[match(c(newnode, node),names(delta))]=1-delta[match(c(newnode, node),names(delta))]
		
		# update values
		new.values=values
		if(direction=="up") {
			new.values=assigndescendants(values, newnode, val, phy, shifts)
			h = ((1/2)*(1/length(get.desc.of.node(newnode,phy)))) / (1/2)
			lnh=log(h)
		} else if(direction=="down") {
			anc=get.ancestor.of.node(node, phy)
			if(anc==root){
				tmp=get.desc.of.node(root, phy)
				dd=tmp[!tmp%in%c(node,shifts)]
				anc.val=unique(as.numeric(values[match(dd, names(values))]))
			} else {
				anc.val=as.numeric(values[which(names(values)==anc)])
			}
			new.values=assigndescendants(values, anc, anc.val, phy, shifts)
			new.values=assigndescendants(new.values, newnode, val, phy, shifts)
			h = (1/2) / ((1/2)*(1/length(get.desc.of.node(node,phy))))
			lnh=log(h)
		} else if(direction=="root"){
			lnh=0
		}
		
#		print(direction)
#		plot(phy, edge.color=ifelse(new.values==1, 1, ifelse(new.values==2, "red", "green")))
#		nodelabels(frame="cir",cex=0.75,bg=NA)
		return(list(new.delta=new.delta, new.values=new.values, lnHastingsRatio=lnh))
	}
}



#proposal utility: move jump-point in tree 
#author: JM EASTMAN 2011

adjustjump <-
function(jump.table, jVCV, k, r, add=FALSE, drop=FALSE, swap=FALSE, phy, node.des){
	
	if(add==TRUE){
		
		jj=tree.slide(phy)
		jump.table=c(jump.table, jj)
		jVCV=updateJvcv(jVCV, from=NULL, to=jj, node.des)
		todo="add"
	} else if(drop==TRUE | swap==TRUE) {
		
		possible.locs=seq(1,length(jump.table))
		adjustedjump.index<-aji<-possible.locs[ceiling(runif(1)*length(possible.locs))]
		adjustedjump.loc<-ajl<-jump.table[adjustedjump.index]
		
		if(drop==TRUE) {
			
			jump.table=jump.table[-adjustedjump.index]
			jVCV=updateJvcv(jVCV, from=adjustedjump.loc, to=NULL, node.des)
			todo="drop"
		} else if(swap==TRUE){
			
			adjustedjump.newloc<-ajnl<-tree.slide(phy, adjustedjump.loc, up=NA)
			jump.table[adjustedjump.index]=adjustedjump.newloc		
			if((from<-adjustedjump.loc) != (to<-adjustedjump.newloc)) {
				jVCV=updateJvcv(jVCV, from=from, to=to, node.des)
			} else {
				jVCV=jVCV
			}
			todo="swap"
		} else {
			stop("Error encountered.")
		}
	} else {
		stop("Must 'add', 'drop', or 'swap' jumps.")
	}
	
	return(list(jump.table=jump.table, jVCV=jVCV))
}



#general phylogenetic utility: recurse down tree changing values of descendants to 'value' until an 'excluded' descendant subtree is reached
#author: JM EASTMAN 2011

assigndescendants <-
function(vv, node, value, phy, exclude=NULL){
	try(vv[match(node, phy$edge[,2])]<-value,TRUE); if(inherits(value, "try-error")) print(list(vv,node,value,exclude))
	desc=get.desc.of.node(node, phy)
	desc=desc[!desc%in%exclude]
	if(length(desc)){
		for(d in desc) vv=assigndescendants(vv, d, value, phy, exclude=exclude)
	}
	return(vv)
}



#general phylogenetic utility for selecting one element of a list of values, each of which can be associated with the edges of a phylogeny (phy$edge[,2])
#author: JM EASTMAN 2011

choose.one <-
function(cur.delta, phy, internal.only=FALSE, edge.prob=0)
# updated 02.26.2011 to weight by branch length (rather than equal weight for each edge)
{
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
	
	# ensure that at least one descendant of root is unshifted 
	root.des=get.desc.of.node(Ntip(phy)+1, phy)
	root.delta=as.logical(bb[match(root.des, phy$edge[,2])])
	if(all(root.delta)) {
		return(choose.one(cur.delta, phy, internal.only, edge.prob)) 
	} else {
		return(bb)
	}
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
function(rates, prop.width, min=0, max=Inf, tuner=0.5, prior.mean) {
	ss=sample(rates, 1)
	ww=which(rates==ss)
	
	if(runif(1)<tuner){
		nn=proposal.slidingwindow(ss, prop.width, min=min, max=max)
	} else {
		nn=proposal.multiplier(ss, prop.width)
	}
	
	nv=nn$v
	lhr=nn$lnHastingsRatio
	lpr=priorratio.exp(ss, 1, nv, 1, prior.mean)

	rates[ww]=nv
	
	return(list(values=rates, lnHastingsRatio=lhr, lnPriorRatio=lpr))

}


#rjmcmc proposal mechanism for updating a single numeric class from a vector
#author: JM EASTMAN 2010

tune.value <-
function(values, prop.width, min=-Inf, max=Inf, tuner=0.5) {
	ss=sample(values, 1)
	ww=which(values==ss)	
	
	if(runif(1)<tuner){
		nn=proposal.slidingwindow(ss, prop.width, min=min, max=max)
	} else {
		nn=proposal.multiplier(ss, prop.width)
	}
	
	nv=nn$v
	lhr=nn$lnHastingsRatio
	lpr=0
	
	values[ww]=nv

	return(list(values=values, lnHastingsRatio=lhr, lnPriorRatio=lpr))
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




#scales VCV matrix of phylogeny by relative evolutionary rates under Brownian motion
#author: JM EASTMAN 2010

updatevcv <-
function(ape.tre, new.rates) {
	n=ape.tre
	n$edge.length=n$edge.length*new.rates
	vv=vmat(n)
	return(vv)
}



updateJvcv=function(jVCV, from=NULL, to=NULL, node.des){
# from: a node for which a jump has been removed
# to: a node for which a jump has been added
	if(!is.null(from)) {
		ff=node.des[[which(names(node.des)==from)]]
		jVCV[ff,ff]=jVCV[ff,ff]-1
	}
	
	if(!is.null(to)){
		tt=node.des[[which(names(node.des)==to)]]
		jVCV[tt,tt]=jVCV[tt,tt]+1
	}
	jVCV
}



updateLvcv=function(Tvcv, Lvcv, sigsqB, sigsqJ){
	return(sigsqB*Tvcv+sigsqJ*Lvcv)
}



jumps.edgewise=function(phy, jump.table){
	jj=rep(0,nrow(phy$edge))
	if(length(jump.table)) {
		jj.tmp=table(jump.table)
		jj[match(names(jj.tmp),phy$edge[,2])]=jj.tmp		
	}
	jj
}




