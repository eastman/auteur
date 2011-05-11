\name{shifts.plot}
\alias{shifts.plot}
\title{plotting of Bayesian posterior samples}
\description{Plots a phylogeny with topological indication of posterior support for evolutionary shifts}
\usage{shifts.plot(phy, base.dir, burnin = 0, level = 0.01, internal.only = FALSE, paint.branches = TRUE, legend=TRUE, pdf = TRUE, verbose = TRUE, lab = "OUT", ...)}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{base.dir}{a directory name, in which the results from \code{rjMCMC} are to be found}
  \item{burnin}{proportion of the chain to be treated as burnin (e.g., from 0 to 1)}
  \item{level}{a level of posterior support (from 0 to 1) for a branch-associated evolutionary shifts, above which results will be compiled}
  \item{internal.only}{dictates whether tips are excluded from processing of evolutionary shifts (\code{internal.only=TRUE})}
  \item{paint.branches}{a logical switch governing whether branches in the phylogeny are to be shaded according to their model-averaged estimates}
  \item{legend}{a logical switch that determines whether keys for the generated symbols and colors are plotted}
  \item{pdf}{a logical switch that determines whether a .pdf of the plot is written to the 'results' subdirectory of the stored Bayesian output}
  \item{verbose}{a logical switch that determines whether additional diagnositics are generated for all branches receiving posterior support beyond 
the defined \code{level}; if \code{verbose=TRUE}, results will be stored in a 'results' subdirectory}
  \item{lab}{a string used to uniquely identify the .pdf if created}
  \item{\dots}{additional arguments to be passed to \code{\link[ape]{plot.phylo}}}
}
\value{
Returns, as a side-effect, a plotted phylogenetic tree with colorized summary of posterior estimates.  Colors that strongly differ from gray are lower (blue) or higher (red) 
than the median posterior estimate across the tree.  Bubbles at branches are proportional to posterior support for a shift occurring at the indicated branch, conditioned on the underlying 
evolutionary process involving at least one shift.}

\author{JM Eastman}

\examples{

# generate tree
n=24
phy=prunelastsplit(birthdeath.tree(b=1,d=0,taxa.stop=n+1))

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
plot(phy,edge.col=ifelse(e==1,cols[1],cols[2]), edge.width=2)
mtext("expected pattern of rates")

## simulate data on the 'rate-shifted' tree
dat=rTraitCont(phy=rphy, model="BM", sigma=sqrt(0.1))

## run two short reversible-jump Markov chains
r=paste(sample(letters,9,replace=TRUE),collapse="")
lapply(1:2, function(x) rjmcmc.bm(phy=phy, dat=dat, ngen=10000, sample.freq=10, prob.mergesplit=0.1, simplestart=FALSE, prop.width=1, fileBase=paste(r,x,sep=".")))

 # collect directories
dirs=sapply(1:2, function(x) dir("./",pattern=paste("BM",r,x,sep=".")))
pool.rjmcmcsamples(base.dirs=dirs, lab=r)

## plot Markov sampled rates
dev.new()
shifts.plot(phy=phy, base=paste(r,"combined.rjMCMC",sep="."), burnin=0.5, legend=TRUE, pdf=FALSE, edge.width=2)

## PASTE UNCOMMENTED FOLLOWING LINE TO DROP DIRECTORIES CREATED BY RJMCMC
 # unlink(dir(pattern=paste(r)),recursive=TRUE)
}
