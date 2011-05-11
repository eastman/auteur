\name{intercalate.samples}
\alias{intercalate.samples}
\title{pool data}
\description{
Combines datasets respecting indices of each}
\usage{intercalate.samples(list.obj)}
\arguments{
  \item{list.obj}{a list of vectors or dataframes}
}
\value{A dataframe of pooled values, whose elements are ordered as given by the order of objects in \code{list.obj}}
\author{JM Eastman}
\examples{
# with character vectors
a=letters[1:10]
A=toupper(letters[1:10])
intercalate.samples(list(a,A))

# with multi-column dataframes
list.obj=lapply(1:5, function(x) {y=data.frame(matrix(data=runif(9)+x, nrow=3)); names(y)=letters[1:3]; return(y)})
intercalate.samples(list.obj)
}