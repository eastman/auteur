\name{calibrate.proposalwidth}
\alias{calibrate.proposalwidth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
initialize proposal width
}
\description{
Generates a reasonable proposal width to initiate sampling for Markov sampling
}
\usage{
calibrate.proposalwidth(phy, dat, nsteps = 100, model, widths = NULL, lim = list(min = 0, max = 1000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
a phylogenetic tree of class 'phylo'}
  \item{dat}{
a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{nsteps}{
number of proposal steps over which to assess proposal widths}
  \item{model}{
currently either \code{model="BM"} or \code{model="OU"} for Brownian motion or an Ornstein-Uhlenbeck process; 
see \code{\link[auteur]{rjmcmc.bm}} for further information on fitting these models}
  \item{widths}{
if unspecified, a series of proposal widths from 1/8 to 8 will be considered}
  \item{lim}{a numeric list defining lower (\code{min}) and upper (\code{max}) bounds on relative rates}

}
\details{
This is primarily an internal function, although may be useful for constraining subsequent runs after an adequate proposal width has been approximated by this function. This 
function is called internally by \code{rjmcmc}, which causes additional sampling (one-tenth total Markov chain length) to be prepended to a run. The sole purpose of this initial 
sampling period is to optimize the proposal width (\code{prop.width}) to ensure that the Markov chain achieves convergence.  Estimates from this calibration are not stored and do not become 
available to the user.}
\author{
JM Eastman}
\examples{
n=4
phy=rescaleTree(phy=rcoal(n=n),totalDepth=100)
dat=rTraitCont(phy=phy, model="BM", sigma=sqrt(0.1))
r=paste(sample(letters,9,replace=TRUE),collapse="")

## with calibrated proposal width
rjmcmc.bm(phy=phy, dat=dat, ngen=5000, sample.freq=10, prob.mergesplit=0.05, prop.width=NULL, fileBase=r)

## with enforced (and large) proposal width
rjmcmc.bm(phy=phy, dat=dat, ngen=5000, sample.freq=10, prob.mergesplit=0.05, prop.width=10, fileBase=r)

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)
}
