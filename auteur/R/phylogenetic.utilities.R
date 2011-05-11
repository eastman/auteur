oumat<-function (vmat, alpha, sigma) 
{
	tips=unique(dim(vmat))
	
	out <- .Call("oumat",	
				 TIMES=list(vmat=as.double(array(vmat))),
				 ALPHA=as.double(alpha^2),
				 SIGMA=as.double(sigma), 
				 TIPS=as.integer(tips),
				 PACKAGE = "auteur")
	return(matrix(out$oumat,ncol=tips))
}



ouweights <- function(phy, epochs, object, ancestors, alpha){
	if(!all(c("anc","des","time","regimes")%in%names(object)))stop("Cannot decipher supplied 'object'.")
	ntip = Ntip(phy)
	object=object[order(object$des),]
	optima = sort(unique(as.numeric(object$regimes)))
    nreg = length(optima)
	
	out <- .Call("ouweights",	
				 PHYLO=list(
							epochs=epochs,
							ancestors=ancestors,
							regimes=as.double(object$regimes),
							optima=as.double(optima),
							alpha=as.double(alpha^2),
							tips=as.integer(ntip)),
				 PACKAGE = "auteur")
	
	W=matrix(out$weights,ncol=nreg)
	opt=optima
	return(list(W=W, optima=opt))
} 

vmat <- function(phy){
	n=Ntip(phy)
	out <- .Call("vmat", tree=list(
									ROOT = as.integer(n+1),
									MAXNODE = as.integer(max(phy$edge[,1])),
									ENDOFCLADE = as.integer(dim(phy$edge)[1]),
									ANC = as.integer(phy$edge[,1]),
									DES = as.integer(phy$edge[,2]),
									EDGES = as.double(c(phy$edge.length,0)),
									VCV = as.double(array(matrix(0, n, n)))),
									PACKAGE = "auteur")
	v=matrix(out$VCV,nrow=n,byrow=FALSE)
	rownames(v)<-colnames(v)<-phy$tip.label
	return(v)
} 

		

ultrametricize=function(phy, tol=1e-8, trim=c("min","max","mean")){
	paths=diag(vmat(phy))
	trim=switch(match.arg(trim),
				min = min(paths),
				max = max(paths),
				mean = mean(paths))
	if(diff(range(paths))<=tol) {
		for(i in 1:length(paths)){
			e=phy$edge.length[which(phy$edge[,2]==i)->bl]
			phy$edge.length[bl]=e+(trim-paths[i])
		}
		return(phy)
	} else {
		stop("Difference in path lengths is larger than supplied tolerance")
	}
}

#general phylogenetic utility for finding the edge labels (from phy$edge[,2]) associated with the most exclusive monophyletic grouping containing at least supplied 'tips' where 'tips' are elements of the phy$tip.label vector
#author: JM EASTMAN 2010

collect.nodes <-
function(tips,phy,strict=FALSE){
	if(length(tips)==1) stop("Must supply at least two tips.")
	tt=phy$tip.label
	tt=lapply(tips, function(x) which(phy$tip.label==x))
	nodes=lapply(tt, function(x) get.ancestors.of.node(x, phy))
	all.nodes=nodes[[1]]
	for(i in 2:length(nodes)) {
		all.nodes=intersect(all.nodes,nodes[[i]])
	}
	base.node=max(all.nodes)
	if(strict) {
		u=unlist(nodes)
		return(c(sort(unique(u[u>=base.node])), unlist(tt)))
	} else {
		return(c(base.node,get.descendants.of.node(base.node,phy)))
	}
}

#general phylogenetic utility for finding the phy$tip.labels associated with the most exclusive monophyletic grouping containing at least supplied 'tips' where 'tips' are elements of the phy$tip.label vector
#author: JM EASTMAN 2010

find.relatives <-
function(tips,phy){
	nn=collect.nodes(tips, phy, strict=FALSE)
	return(phy$tip.label[nn[nn<=Ntip(phy)]])
}

#general phylogenetic utility for returning first ancestor (as a numeric referent) of the supplied node 
#author: JM EASTMAN 2010

get.ancestor.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,2]==node),1])
}



#general phylogenetic utility for returning all ancestors (listed as given in phy$edge[,2]) of a node 
#author: JM EASTMAN 2010

get.ancestors.of.node <-
function(node, phy) {
	a=c()
	if(node==(Ntip(phy)+1->root)) return(NULL)
	f=get.ancestor.of.node(node, phy)
	a=c(a,f)
	if(f>root) a=c(a, get.ancestors.of.node(f, phy))
	return(a)
}


#general phylogenetic utility for returning the first (usually, unless a polytomy exists) two descendants the supplied node 
#author: JM EASTMAN 2010

get.desc.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}


#general phylogenetic utility for returning all descendants (listed as given in phy$edge[,2]) of a node (excluding the supplied node) 
#author: JM EASTMAN 2010, based on GEIGER:::node.leaves()

get.descendants.of.node <-
function (node, phy, tips=FALSE) 
{
    node <- as.numeric(node)
    n <- Ntip(phy)
	if(node <= n) return(NULL)
    l <- c()
    d <- get.desc.of.node(node, phy)
    for (j in d) {
		l <- c(l, j)
        if (j > n) {
			l <- c(l, get.descendants.of.node(j, phy))
		}
    }
    if(tips) return(l[l<=n]) else return(l)
}


get.epochs<-function(phy) 
# ordered: tip to root
{
	ancestors=lapply(1:Ntip(phy), function(x) c(x,rev(sort(get.ancestors.of.node(x,phy)))))
	phylo=get.times(phy)
	tt=phylo$time[order(phylo$des)]
	e=lapply(ancestors, function(x) tt[x])
	names(e)=phy$tip.label
	return(e)
}

get.times<-function(phy) 
{
	df=data.frame(phy$edge)
	names(df)=c("anc","des")
	anc=lapply(df$des,get.ancestors.of.node,phy)
	order.edges=c(0,phy$edge.length)[order(c(Ntip(phy)+1,df$des))]
	df$time=0
	for(n in 1:length(anc)){
		df$time[n]=sum(order.edges[anc[[n]]])+order.edges[df$des[n]]
	}
	
	df=rbind(c(0,Ntip(phy)+1,0),df)
	
	return(df)
}

#general phylogenetic utility for determining whether a node is the root of the phylogeny
#author: JM EASTMAN 2010

is.root <-
function(node,phy) {
	if(node==phy$edge[1,1]) return(TRUE) else return(FALSE)
}

#general phylogenetic utility for returning all ancestors of nodes in a phylogeny
#author: JM EASTMAN 2011

get.ancestors <-
function(phy) {
	r=lapply(1:Ntip(phy), function(x) c(x,rev(sort(get.ancestors.of.node(x,phy)))))
	r
}


#general phylogenetic utility for pruning the most recent speciation event from a phylogeny
#author: JM EASTMAN 2010

prunelastsplit<-function(phy) {
	require(ape)
	nn=as.numeric(names(xx<-branching.times(phy))[which(xx==min(xx))])
	for(node in 1:length(nn)) {	
		sp=sample(phy$edge[which(phy$edge[,1]==nn[node]),2],1)
		obj=drop.tip(phy,phy$tip.label[sp])
	}
	return(obj)
}

