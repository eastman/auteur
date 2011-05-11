#rjmcmc proposal mechanism: proposal mechanism for traversing model complexity (by one parameter at a time)
#author: JM EASTMAN 2010

splitormerge <-
function(cur.delta, cur.values, phy, node.des, lambda=lambda, theta=FALSE, internal.only=FALSE) { 
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
		if(theta) u=splitvalue(cur.vv=cur.vv, n.desc=n.desc, n.split=n.split, factor=lambda) else u=splitrate(cur.vv=cur.vv, n.desc=n.desc, n.split=n.split, factor=lambda)
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
		lnPriorRatio = log(ptpois(K+1,lambda,nn)/ptpois(K,lambda,nn))
		
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
		lnPriorRatio = log(ptpois(K-1,lambda,nn)/ptpois(K,lambda,nn))
		
	}
	
	return(list(new.delta=new.bb, new.values=new.vv, lnHastingsRatio=lnHastingsRatio, lnPriorRatio=lnPriorRatio, decision=decision))
}

