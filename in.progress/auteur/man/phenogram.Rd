\name{phenogram}
\alias{bm.phenogram}
\alias{levy.phenogram}
\alias{phenogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
trait evolution along phylogenies}
\description{
simulating incremental phenotypic evolution by Brownian motion or Levy process along a phylogenetic tree}
\usage{
phenogram(phy, model = c("levy", "bm"), alpha = 0, rate = 0.01, lambda = 0.5, sigma = 0.5, increments = 100, plot = TRUE, node.cex = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
a phylogenetic tree of class \code{phylo}}
  \item{model}{
a string (either "levy" or "bm") determining which process to use in simulating trait evolution}
  \item{alpha}{
the ancestral phenotypic value}
  \item{rate}{
a numeric value determining the rate of evolution; rates are not comparable between Levy and Brownian motion}
  \item{lambda}{
mean of the Poisson distribution from which the number of trait-evolutionary saltations are drawn; applies only where \code{model="levy"}}
  \item{sigma}{
standard deviation of the distribution from which saltation sizes are drawn; applies only where \code{model="levy"}}
  \item{increments}{
the number of branch segments used to discretize the time-continuous process (for a root-to-tip path); if \code{increments=NULL}, only 
the phenotypic values at nodes are plotted}
  \item{plot}{a logical value determining whether a plot is generated}
  \item{node.cex}{plot symbol size for the phenotypic values at nodes}
  \item{\dots}{
additional arguments to be passed to \code{plot}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function generates a phenogram (if \code{plot=TRUE}) of evolutionary change in phenotype along a phylogeny and returns a table of 
ancestral values (\code{phenotype}) for each node in the topology (\code{descendant}). Using the \code{\link[ape]{ape-package}} format 
for edge identification, each \code{ancestor} and \code{descendant} is given an integer value. The split \code{time} of each node 
is given in this last column of the returned table.  }
\references{
%% ~put references to the literature/web site here ~
}
\author{
LJ Revell, LJ Harmon, and JM Eastman}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## BROWNIAN MOTION ##
phenogram(rcoal(20)->phy,model="bm",rate=0.005,increments=100,lwd=2)
phenogram(phy,model="bm",rate=0.0075,increments=100,lwd=2)
phenogram(phy,model="bm",rate=0.01,increments=100,lwd=2)
phenogram(phy,model="bm",rate=0.01,increments=NULL,lwd=2)

## LEVY MOTION ##
phenogram(phy,model="levy",rate=0.1,lambda=0.25,sigma=0.75,increments=100,lwd=2)
phenogram(phy,model="levy",rate=0.1,lambda=0.5,sigma=0.75,increments=100,lwd=2)
phenogram(phy,model="levy",rate=0.1,lambda=2,sigma=0.75,increments=100,lwd=2)
phenogram(phy,model="levy",rate=0.1,lambda=2,sigma=0.75,increments=NULL,lwd=2)

}