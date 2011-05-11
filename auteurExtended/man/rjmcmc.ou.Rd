\name{rjmcmc.ou}
\alias{rjmcmc.ou}
\title{Bayesian sampling of shifts in trait optima: Ornstein-Uhlenbeck process}
\description{Implements reversible-jump Markov chain Monte Carlo sampling to identify shifts in the process of continuous trait evolution across phylogenetic trees}
\usage{rjmcmc.ou(phy, dat, SE = 0, ngen = 1000, sample.freq = 100, prob.mergesplit = 0.1, prob.theta = 0.4, prob.alpha = 0.2, upper.alpha = 100, lambdaK = log(2), constrainK = FALSE, prop.width = NULL, simplestart = FALSE, internal.only = FALSE, summary = TRUE, fileBase = "result")}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}; see \code{\link[geiger]{name.check}}}
  \item{SE}{a named vector of standard errors for each trait value; applied to all trait values if given a single value}
  \item{ngen}{number of sampling generations}
  \item{sample.freq}{frequency with which Markov samples are retained (e.g., \code{sample.freq=10} retains every tenth sample in the chain)}
  \item{prob.mergesplit}{proportion of proposals that split or merge optima, thereby altering complexity of the evolutionary model}
  \item{prob.theta}{proportion of proposals in which selective regimes for branches are updated}
  \item{prob.alpha}{proportion of proposals in which the evolutionary constraint parameter is updated}
  \item{upper.alpha}{an upper bound on the evolutionary constraint parameter}
  \item{lambdaK}{the shape parameter for the Poisson prior-distribution on number of distinct optima}
  \item{constrainK}{a constraint of model complexity (e.g., \code{constrainK=2} restricts sampling to models with two independent optima)}
  \item{prop.width}{if \code{NULL}, the proposal width is calibrated; otherwise, proposal width is constrained to the supplied value}
  \item{simplestart}{if \code{FALSE}, a random model complexity is chosen; if \code{TRUE}, Markov sampling begins under a global-optimum model}
  \item{internal.only}{dictates whether any branch in the tree can be chosen for a regime shift or tips only (\code{internal.only=FALSE})} 
  \item{summary}{a logical switch that determines whether any output is printed to the console after the analysis has finished}
  \item{fileBase}{a string used to uniquely identify the base directory to which Markov samples are written}
}
\details{
In the two primary objects that store posterior estimates of \code{shifts} and \code{optima} (see \bold{Value}), branches are labeled in accordance with the numeric identifiers found in \code{phy$edge}, which represents a table of ancestor-descendant
 relationship (see \code{\link[ape]{read.tree}}.  In the \bold{Examples} below do \bold{not} expect convergence with such short chains!
}
\value{
If \code{summary=TRUE} the global rate of acceptance and calibrated proposal width is returned as a list.  After a run has completed, the rates of acceptance for various proposal mechanisms are printed, along with 
the (calibrated) proposal width used for chain sampling. 

Posterior results are written to several files within a base directory, the contents of which are as follows:
  \item{logfile}{a logfile including the following for each Markov sample; the generations at which samples were retained; the generalized model of trait evolution used for inference (e.g., Brownian motion or 'BM"); 
  the mean treewide optimum trait value in the sampled generation; the number of independent optima; and the likelihood of the model at each sample.  Do not be alarmed if it seems to take some time before values 
  begin to appear in the \code{logfile}: samples are not saved until after the calibration period has terminated }
  \item{errors}{a logfile to which errors are recorded, if any are generated in the Markov chain}
  \item{posteriorsamples}{a compressed .rda file, storing estimates from the Markov chain for relevant parameters; contents of this file can be loaded 
  with \code{load("posteriorsamples.rda")}. This list object, whose referent is also \code{posteriorsamples} when in R-memory, includes two dataframes: 
  \code{optima} comprises the inferred selective regime (trait optimum) for each branch in the phylogeny; \code{shifts} stores the branches that were chosen for 
 shifts in selective regime for each retained sample of the Markov chain (\code{1} signifies a shift; 
  \code{0} implies a branch that was not sampled for a regime shift). Column names of both components of \code{posteriorsamples} are from the numeric vector of 
  branch labels (see \bold{Details}). For further information on .Rd and .rda files, see \code{\link[base]{save}} and \code{\link[base]{load}} as well as the 
  \bold{Examples} below. }
}
\author{AA King and JM Eastman}
\note{This function is largely based on code from \code{\link[ouch]{ouch-package}}}

\examples{
#############
## generate tree
n=24
while(1) {
	phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))
	phy$tip.label=paste("sp",1:n,sep="")
	rphy=reorder(phy,"pruningwise")

	# find an internal edge
	anc=get.desc.of.node(Ntip(phy)+1,phy)
	branches=phy$edge[,2]
	branches=branches[branches>Ntip(phy) & branches!=anc]
	branch=branches[sample(1:length(branches),1)]
	desc=get.descendants.of.node(branch,phy)
	if(length(desc)>=4) break()
}
rphy=phy

## build array of optima
t=numeric(nrow(phy$edge))
t[match(c(branch,desc),phy$edge[,2])]=8
cols=c("red","gray")
dev.new()
plot(phy,edge.col=ifelse(t==8,cols[1],cols[2]), edge.width=2) 
mtext("expected pattern of optima")

#############
## simulate data on the 'optimum-shifted' tree: shift from optimum 0.0 to optimum 2.0
dat=rTraitCont(phy=rphy, model="OU", sigma=sqrt(0.01), alpha=0.1, theta=t)

## run two short reversible-jump Markov chains
r=paste(sample(letters,9,replace=TRUE),collapse="")
lapply(1:2, function(x) rjmcmc.ou(phy=phy, dat=dat, ngen=10000, sample.freq=10, prob.mergesplit=0.1, prob.theta=0.4, prob.alpha=0.2, simplestart=FALSE, prop.width=NULL, fileBase=paste(r,x,sep=".")))

 # collect directories
dirs=dir("./",pattern=paste("OU",r,sep="."))
pool.rjmcmcsamples(base.dirs=dirs, lab=r)

## plot Markov sampled optima
dev.new()
shifts.plot(phy=phy, base=paste(r,"combined.rjmcmc",sep="."), burnin=0.5, legend=TRUE, pdf=FALSE, edge.width=2) 

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)
 }