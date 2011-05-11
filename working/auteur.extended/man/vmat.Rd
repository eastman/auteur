\name{vmat}
\alias{vmat}
\title{computation of phylogenetic variance-covariance matrix}
\description{Calculates the VCV matrix for a phylogenetic tree}
\usage{vmat(phy)}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
}
\details{
This function is a conversion of \code{\link[ape]{vcv.phylo}} into compiled \code{C++} for rapid generation of the expected trait-variances 
and trait covariances among species under Brownian motion evolution.  This function is efficient for trees with fewer than ca. 5000 tips; for larger trees, 
sufficient memory (> 2 Gb RAM) may be required.  
}
\value{
A variance-covariance matrix for all tips within the supplied tree.
}
\author{JM Eastman, based on \code{\link[ape]{vcv.phylo}} by Emmanuel Paradis}
\examples{
## generate tree
n=250
phy=rescaleTree(phy=rcoal(n=n),totalDepth=100)

## compare function times for vcv.phylo() and vmat()
print(system.time(vcv.phylo(phy)))
print(system.time(vmat(phy)))

## generate some smaller matrices
n=4
phy=rescaleTree(phy=rcoal(n=n),totalDepth=100)

## compute the variance-covariance matrix with ape and rjmcmc
vcv.phylo(phy)->vAPE
vmat(phy)->vRJ

## print the matrices
print(vAPE)
print(vRJ)

## verify that both packages return identical results
print(all(vAPE==vRJ))
}
