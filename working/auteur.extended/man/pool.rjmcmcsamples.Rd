\name{pool.rjmcmcsamples}
\alias{pool.rjmcmcsamples}
\title{combining posterior samples from Bayesian analysis}
\description{Generates a pooled sample of posterior estimates from any number of independent runs}
\usage{pool.rjmcmcsamples(base.dirs, lab = "")}
\arguments{
  \item{base.dirs}{
a vector of directory names, in which the results from \code{rjmcmc} are to be found}
  \item{lab}{
an optional string used to name the directory to which pooled results are written}
}

\examples{
n=24
phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))
dat=rTraitCont(phy=phy, model="BM", sigma=sqrt(0.1))

## run three short reversible-jump Markov chains
r=paste(sample(letters,9,replace=TRUE),collapse="")
lapply(1:3, function(x) rjmcmc.bm(phy=phy, dat=dat, ngen=500, sample.freq=10, prob.mergesplit=0.1, simplestart=TRUE, jumpsize=1, fileBase=paste(r,x,sep=".")))

 # collect directories
dirs=sapply(1:3, function(x) dir("./",pattern=paste("BM",r,x,sep=".")))
pool.rjmcmcsamples(base.dirs=dirs, lab=r)

 # plot density of sampled mean rates across tree
rates=read.table(file=paste(dir(pattern=paste(r,"combined.rjmcmc",sep=".")),paste(r,"rjmcmc.log",sep="."),sep="/"),header=TRUE)$mean.rate
trace.plot(rates,col="maroon",alpha=0.4,xlab="rates",log="x",legend.control=list(plot=FALSE))

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)
}
