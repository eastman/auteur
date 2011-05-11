\name{calibrate.jumpsize}
\alias{calibrate.jumpsize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
initialize proposal width
}
\description{
Generates a reasonable proposal width to initiate sampling for Markov sampling
}
\usage{
calibrate.jumpsize(phy, dat, nsteps = 100, jumpsizes = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
a phylogenetic tree of class 'phylo'}
  \item{dat}{
a named vector of continuous trait values, associated with each species in \code{phy}}
  \item{nsteps}{
number of proposal steps over which to assess proposal widths}
  \item{jumpsizes}{
if unspecified, a series of jumpsizes from 1/8 to 8 will be considered}
}
\details{
This is primarily an internal function, although may be useful for constraining subsequent runs after an adequate proposal width has been approximated by this function.}
\author{
JM Eastman}
\examples{
n=4
phy=rescaleTree(phy=rcoal(n=n),totalDepth=100)
dat=rTraitCont(phy=phy, model="BM", sigma=sqrt(0.1))
r=paste(sample(letters,9,replace=TRUE),collapse="")

## with calibrated jumpsize
rjmcmc.bm(phy=phy, dat=dat, ngen=5000, sample.freq=10, prob.mergesplit=0.05, jumpsize=NULL, fileBase=r)

## with enforced (and large) jumpsize
rjmcmc.bm(phy=phy, dat=dat, ngen=5000, sample.freq=10, prob.mergesplit=0.05, jumpsize=10, fileBase=r)

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)
}
