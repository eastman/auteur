\name{trace.plot}
\alias{trace.plot}
\title{plotting of Bayesian posterior samples}
\description{Compares posterior densities of Bayesian estimates}
\usage{trace.plot(obj, col, alpha, lwd = 1, hpd = 0.95, bars = TRUE, legend.control = list(plot = TRUE, pos = NA, cex = 1, pt.cex = 1, pch = 22, title = ""), truncate = list(min = NA, max = NA), xlim = list(min = NA, max = NA), ...)}
\arguments{
  \item{obj}{a vector or dataframe (where each column is separately plotted); 
  names of columns are retrieved if a legend is to be generated}
  \item{col}{either empty or a vector of colors with one for each column in the \code{obj}}
  \item{alpha}{a degree of color translucence to be used for shaded densities, where 1 is opaque, 0 is clear}
  \item{lwd}{line width used for highest-density range and density outlines (see \code{\link{par}})}
  \item{hpd}{either \code{NULL} or a value between 0 and 1, corresponding to the width of the highest 
  density region to be shaded}
  \item{bars}{a logical specifier which, if \code{bars=TRUE}, outlines the width of the high density region}
  \item{legend.control}{if \code{plot=TRUE}, a legend is generated (see \code{\link{legend}} for details concerning this list)}
  \item{truncate}{if \code{min} and (or) \code{max} are defined, values below and (or) above the truncation 
  point are removed prior to computation of the highest density region and prior to plotting}
  \item{xlim}{if \code{min} and (or) \code{max} are defined, values below and (or) above the limit(s) are 
  excluded from the plot without further corrupting the data}
  \item{\dots}{further arguments to be passed to \code{\link{plot}}}
}
\author{JM Eastman}

\seealso{\code{\link[diversitree]{profiles.plot}}, upon which this function is largely based}
\examples{
# construct and plot a dataframe of three 'traces', excluding the largest values
d=data.frame(list(rnorm(n=1000,mean=9,sd=5)),list(rnorm(n=1000,mean=1.5,sd=2)),list(rnorm(n=1000,mean=-0.5,sd=1.5)))
names(d)=letters[1:3]
trace.plot(obj=d, xlim=list(max=10), col=c("maroon","gray","purple"), alpha=0.3, lwd=2)

dev.new()
# construct and plot a dataframe of two 'traces' using a log scale
d=data.frame(list(rlnorm(n=1000,meanlog=1,sdlog=0.5)),list(rlnorm(n=1000,meanlog=1.5,sdlog=0.6)))
names(d)=letters[1:2]
trace.plot(obj=d,col=c("maroon","purple"),alpha=0.3,log="x")
}
