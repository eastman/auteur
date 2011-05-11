\name{compare.rates}
\alias{compare.rates}
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
n=24
phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))

## simulate a rate shift
 # find an internal edge
anc=get.desc.of.node(Ntip(phy)+1,phy)
branches=phy$edge[,2]
branches=branches[branches>Ntip(phy) & branches!=anc]
branch=branches[sample(1:length(branches),1)]
desc=get.descendants.of.node(branch,phy)
rphy=phy
rphy$edge.length[match(desc,phy$edge[,2])]=phy$edge.length[match(desc,phy$edge[,2])]*64
e=numeric(nrow(phy$edge))
e[match(c(branch,desc),phy$edge[,2])]=1
cols=rev(diverge_hcl(n=2))
dev.new()
plot(phy,edge.col=ifelse(e==1,cols[1],cols[2]), edge.width=2)
mtext("expected pattern of rates")

## simulate data on the 'rate-shifted' tree
dat=rTraitCont(phy=rphy, model="BM", sigma=sqrt(0.1))

## run reversible-jump MCMC for a short chain
r=paste(sample(letters,9,replace=TRUE),collapse="")
rjmcmc.bm(phy=phy, dat=dat, ngen=10000, sample.freq=10, prob.mergesplit=0.1, fileBase=r)

 # collect posterior rates 
load(paste(paste("BM",r,"parameters",sep="."),paste(r,"posteriorsamples.rda",sep="."),sep="/"))
rates=posteriorsamples$rates
burnin=round(0.5*nrow(rates))
rates=rates[-c(1:burnin),]

# compare the shifted branches to the background
shift=desc
branches=1:nrow(phy$edge)
background=branches[is.na(match(branches,shift))]

# exclude stem branch from background
background=background[-(which(background==(n+1)))]
comp=compare.rates(branches=list(A=shift, B=background), phy=phy, posterior.values=rates)

# plot posterior distributions of scaled rates by group
prates=data.frame(comp$scl.A, comp$scl.B)
names(prates)=c("shift","background")
trace.plot(obj=data.frame(comp$scl.A, comp$scl.B), col=c("red","black"), alpha=0.4, log="x")

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
# unlink(dir(pattern=paste(r)),recursive=TRUE)
}
