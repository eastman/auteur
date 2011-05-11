\name{randomization.test}
\alias{randomization.test}
\title{statistical comparison of sets of values by randomization}
\description{Compares means by bootstrap resampling of differences between empirical distributions}
\usage{
randomization.test(obs = obs, exp = exp, mu = 0, iter = 10000, two.tailed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{a vector of numeric values}
  \item{exp}{a vector of numeric values}
  \item{mu}{the true difference in means}
  \item{iter}{number of randomization comparisons to perform}
  \item{two.tailed}{as default, the test is performed under a one-tailed assumption; if \code{two.tailed=FALSE}, probability values associated with either tail of the comparison distribution are returned, 
otherwise, a two-tailed result is returned}
}
\details{
If a single value is supplied for \code{obs}, this test equates to finding the quantile in \code{exp} in which \code{obs} would be found (under a one-tailed test); 
see \bold{Examples} and also \code{\link[stats]{ecdf}}}
\value{
A list, whose contents are determined by the above argument: 
  \item{unnamed value}{if \code{two.tailed=TRUE}, this is the two-tailed p-value}
  \item{diffs}{the full resampling distribution of differences between \code{obs} and \code{exp}, given \code{mu} }
  \item{greater}{if \code{two.tailed=FALSE}, this is the p-value associated with the righthand tail}
  \item{lesser}{if \code{two.tailed=FALSE}, this is the p-value associated with the lefthand tail}
}
\author{JM Eastman}
\examples{

# a comparison between two distributions
a=rnorm(n=1000, mean=1, sd=0.5)
b=rnorm(n=1000, mean=0, sd=1)
randomization.test(obs=a, exp=b, two.tailed=FALSE)

# a comparison of a single value to a normal distribution
a=3
b=rnorm(n=1000, mean=0, sd=1)
randomization.test(obs=a, exp=b, two.tailed=FALSE)

# compare above result with ecdf(), in which we compute an empirical 
f=ecdf(b)
print(1-f(a))		# analogous to a one-tailed test as above
}