\name{ou.lik.fn}
\alias{ou.lik.fn}
\title{statistical comparison of posterior rate estimates}
\description{Conducts randomization test to determine if posterior rate estimates for two groups of phylogenetic branches differ}
\usage{compare.rates(branches = list(A = c(NULL, NULL), B = c(NULL, NULL)), phy, posterior.values, ...)}
\arguments{
  \item{branches}{a list of numeric labels for branches in the phylogeny, consistent with the node 
  labels found in \code{phy$edge} (see \code{\link[ape]{read.tree}})}
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{posterior.values}{a dataframe of posterior rate estimates, where rows are sampled 
  generations in the Markov chain and columns distinguish branches in the phylogeny}
  \item{\dots}{additional arguments to be passed to \code{\link{randomization.test}}}}
\details{
For each generation in the Markov sample, a scaled relative rate is computed for each collection of branches (so defined by the supplied list 
of \code{branches}).  Scaling occurs by weighting posterior rate estimates by the length of the each branch in the group.  In so doing, longer 
branches contribute greater weight to the scaled measure of evolutionary rate in each group. Values need not be rates, but elements must be capable of being 
associated with each branch in the phylogeny.}
\value{
A list of scaled relative values for each branch and the result of statistical comparison:
  \item{scl.A }{scaled posterior rate estimates from each Markov sample for group 'A'}
  \item{scl.B }{scaled posterior rate estimates from each Markov sample for group 'B'}
  \item{r.test}{probability value of randomization comparison between groups for phylogenetically-scaled values, as default under a two-tailed test}
}
\author{JM Eastman}
\examples{
## generate tree
	n=400
	phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))
	phy$tip.label=paste("sp",1:n,sep="")
	dat=rTraitCont(phy)
	
	object=get.times(phy)
	r=sample(c(runif(1)),nrow(object),TRUE)
	a=lapply(1:Ntip(phy), function(x) c(x,rev(sort(get.ancestors.of.node(x,phy)))))
	e=get.epochs(phy)
	vcv=vmat(phy)
	rate=runif(1)
	alpha=runif(1)
	system.time(ou.lik.fn(dat, phy, e, object, a, vcv, rate, alpha, r)->lik)
	
	
	## CHECK WITH OUCH ##
	tree=ape2ouch(phy,scale=F)
	dat=dat[order(match(names(dat),tree@nodelabels))]
	otd=as(tree,"data.frame")
	dd=data.frame(labels=names(dat),trait=dat)
	otd=merge(dd, otd, by="labels", all=TRUE)
	otd=otd[otd$nodes,]
	match(seq(1,nrow(otd)),otd$nodes)->m
	otd=otd[m,]
	rownames(otd)=1:nrow(otd)
	
	
	system.time(ouch.ou.lik.fn(tree, alpha^2, as.matrix(rate), otd["trait"], r, unique(r))->oo)

	if(round(oo$lnL-lik$lnL, digits=10)==0) print("MATCH") else print(paste("UNMATCHED: by ", abs(oo$lnL-lik$lnL),sep=""))

	

	
}