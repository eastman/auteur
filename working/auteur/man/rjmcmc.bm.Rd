\name{rjmcmc.bm}
\alias{rjmcmc.bm}
\title{Bayesian sampling of trait evolutionary rates}
\description{Implements reversible-jump Markov chain Monte Carlo sampling to identify shifts in the process of continuous trait evolution across phylogenetic trees}
\usage{rjmcmc.bm(phy, dat, SE=0, ngen = 1000, sample.freq = 100, prob.mergesplit = 0.05, prob.root = 0.05, lambdaK = log(2), constrainK = FALSE, jumpsize = NULL, simplestart = FALSE, internal.only = FALSE, summary = TRUE, fileBase = "result")}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}; see \code{\link[geiger]{name.check}}}
  \item{SE}{a named vector of standard errors associated with each trait value; if a single value is given (e.g., \code{SE=0.05}), the standard error is assumed to apply for all trait values}
  \item{ngen}{number of sampling generations}
  \item{sample.freq}{frequency with which Markov samples are retained (e.g., \code{sample.freq=10} retains every tenth sample in the chain)}
  \item{prob.mergesplit}{a value from 0 to 1, signifying the proportion of proposals that split rates or merge rates, thereby increasing or decreasing complexity of the evolutionary model}
  \item{prob.root}{a value from 0 to 1, governing the proportion of proposals in which the root state is updated}
  \item{lambdaK}{the shape parameter for the Poisson distribution used as a prior on the number of distinct rates in the tree (see \code{\link[stats]{rpois}})}
  \item{constrainK}{if given an integer, Markov sampling becomes dimensionally constrained (e.g., \code{constrainK=2} would restrict posterior 
  samples to be drawn only from models with two independent rates)}
  \item{jumpsize}{if \code{NULL}, the proposal width is intialized through a calibration period (see \code{\link{calibrate.jumpsize}}); otherwise, 
  proposal width is constrained to the supplied value}
  \item{simplestart}{if \code{FALSE}, a random starting point (with \code{k} rate parameters) is chosen; if \code{TRUE}, Markov sampling begins under a global-rate model}
  \item{internal.only}{a logical switch, dictating whether any branch in the tree can be chosen for a rate shift (\code{internal.only=FALSE}) 
  or whether all but tips can be sampled for a rate shift}
  \item{summary}{a logical switch that determines whether any output is printed to the console after the analysis has finished}
  \item{fileBase}{a string used to uniquely identify the base directory to which Markov samples are written}
}

\details{
In the two primary objects that store posterior estimates of \code{rate.shifts} and \code{rates} (see \bold{Value}), branches are labeled in accordance with the numeric identifiers found in \code{phy$edge}, which represents a table of ancestor-descendant
 relationship (see \code{\link[ape]{read.tree}} and the first plot from the \bold{Examples} below).  In the \bold{Examples} below do \bold{not} expect convergence with such short chains!
}
\value{
If \code{summary=TRUE} the global rate of acceptance and calibrated proposal width is returned as a list.  After a run has completed, the rates of acceptance for various proposal mechanisms are printed, along with 
the calibrated jumpsize used for the last nine-tenths of the chain. 

Posterior results are written to several files within a base directory, the contents of which are as follows:
  \item{logfile}{a logfile including the following for each Markov sample; the generations at which samples were retained; the generalized model of trait evolution used for inference (e.g., Brownian motion or 'BM"); 
  the mean treewide rate in the sampled generation; the number of independent rates; the root state; and the likelihood of the model at each sample.  Do not be alarmed if it seems to take some time before values 
  begin to appear in the \code{logfile}: samples are not saved until after the calibration period has terminated }
  \item{errors}{a logfile to which errors are recorded, if any are generated in the Markov chain}
  \item{posteriorsamples}{a compressed .rda file, storing estimates from the Markov chain for relevant parameters; contents of this file can be loaded 
  with \code{load("posteriorsamples.rda")}. This list object, whose referent is also \code{posteriorsamples} when in R-memory, includes two dataframes: 
  \code{rates} comprises the relative-rate estimates for each branch in the phylogeny; \code{rate.shifts} stores the branches that were chosen for 
  rate shifts in each retained sample of the Markov chain (\code{1} signifies a shift; 
  \code{0} implies a branch that was not sampled for a rate shift). Column names of both components of \code{posteriorsamples} are from the numeric vector of 
  branch labels (see \bold{Details}). For further information on .Rd and .rda files, see \code{\link[base]{save}} and \code{\link[base]{load}} as well as the 
  \bold{Examples} below. }
}

\author{LJ Harmon, AL Hipp, and JM Eastman}
\examples{
#############
## generate tree
n=24
phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))

# plot phylogeny with branches labeled as they are identified in by rjmcmc.bm
plot(phy, label.offset = strwidth(par("pch"), cex = 0.35)) ## PLOT ONE ##
edgelabels(text=phy$edge[,2], frame="c", bg=add.transparency("white",0.75), cex=0.75)
mtext("branch labels used by the reversible-jump sampler")

#############
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
plot(phy,edge.col=ifelse(e==1,cols[1],cols[2]), edge.width=2) ## PLOT TWO ##
mtext("expected pattern of rates")

## simulate data on the 'rate-shifted' tree
dat=rTraitCont(phy=rphy, model="BM", sigma=sqrt(0.1))

## run two short reversible-jump Markov chains
r=paste(sample(letters,9,replace=TRUE),collapse="")
lapply(1:2, function(x) rjmcmc.bm(phy=phy, dat=dat, ngen=10000, sample.freq=10, prob.mergesplit=0.1, simplestart=TRUE, jumpsize=1, fileBase=paste(r,x,sep=".")))

 # collect directories
dirs=dir("./",pattern=paste("BM",r,sep="."))
pool.rjmcmcsamples(base.dirs=dirs, lab=r)

## view contents of .rda
load(paste(paste(r,"combined.rjmcmc",sep="."),paste(r,"posteriorsamples.rda",sep="."),sep="/"))
print(head(posteriorsamples$rates))
print(head(posteriorsamples$rate.shifts))

## plot Markov sampled rates
dev.new()
shifts.plot(phy=phy, base=paste(r,"combined.rjmcmc",sep="."), burnin=0.5, legend=TRUE, pdf=FALSE, edge.width=2) ## PLOT THREE ##

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)
}
