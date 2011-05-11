\name{tracer}
\alias{tracer}
\title{plotting of Bayesian posterior samples }
\description{
Generates diagnostics for reversible-jump Markov sampling densities}
\usage{tracer(base.log, lambdaK = log(2), Bayes.factor = NULL, burnin = 0.25, col.line = 1, pdf = TRUE, factor="M",...)}
\arguments{
  \item{base.log}{a path to the stored summary logfile of the Markov sample}
  \item{lambdaK}{the shape parameter for the Poisson distribution used as a prior on the number of distinct rates in the tree (see \code{\link[stats]{rpois}})}
  \item{Bayes.factor}{a vector of two numbers, specifying a Bayes factor comparison between \code{j}- and \code{k}-rate models; alternatively, the argument 
can be \code{Bayes.factor=c(1,"multi")} in which a comparison between a global-rate and multiple-rate model is made}
  \item{burnin}{proportion of the chain to be treated as burnin (e.g., from 0 to 1); burnin is used for all plots but the trace of the log-likelihoods}
  \item{col.line}{a color to be used for density plots (both for a \code{Bayes.factor} comparison and for the density of mean rates across the Markov sample)}
  \item{pdf}{a logical switch that determines whether plots are written to .pdf (and stored within the stored Bayesian output)}
  \item{factor}{either a character ("k" or "M") or \code{NULL}, which determines how likelihood trace axis is defined: in thousands (\code{k}) or millions (\code{M}) of generations}
  \item{\dots}{further arguments to be passed to \code{\link{plot}}}
}
\value{
if \code{pdf=TRUE}, a .pdf will be generated that includes the trace of log-likelihoods, the density of sampled mean rates, the 
correspondence between prior and posterior weights for \code{k}-rate models, and a Bayes factor comparison }
\author{JM Eastman}
\examples{
# generate tree
n=24
phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))

# simulate data 
dat=rTraitCont(phy=phy, model="BM", sigma=sqrt(0.1))

# run reversible-jump MCMC for a short chain
r=paste(sample(letters,9,replace=TRUE),collapse="")
rjmcmc.bm(phy=phy, dat=dat, ngen=10000, sample.freq=10, prob.mergesplit=0.1, fileBase=r)

# plot Markov sampled rates and other traces
tracer(base.log=paste(paste("BM",r,"parameters",sep="."),paste("BM",r,"rjmcmc.log",sep="."),sep="/"), lambdaK=log(2), Bayes.factor=c(1,"multi"), burnin=0.25, pdf=FALSE, factor="k")

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)
}
