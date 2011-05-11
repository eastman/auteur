\name{ouweights}
\alias{ouweights}
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
	require(ouch)

	n=80
	phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))
	phy$tip.label=paste("sp",1:n,sep="")
	
	object=get.times(phy)
	r=sample(c(0,1,-1),nrow(object),TRUE)
	object$regimes=r
	
	a=lapply(1:Ntip(phy), function(x) c(x,rev(sort(get.ancestors.of.node(x,phy)))))
	
	e=get.epochs(phy)
	
	ouweights(phy, e, object, a, alpha<-runif(1))->b; print(b)
	
	
	
}