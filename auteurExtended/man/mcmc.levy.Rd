\name{mcmc.levy}
\alias{mcmc.levy}
\title{Bayesian sampling of shifts in trait-evolutionary rates: jump-diffusion}
\description{Implements Markov chain Monte Carlo sampling to identify shifts in the process of continuous trait evolution across phylogenetic trees}
\usage{mcmc.levy(phy, dat, SE=0, ngen = 1000, sample.freq = 100, prob.lambda = 0.75, prob.sigsqB = 0.075, prob.sigsqJ = 0.15, prop.width = NULL, lambdaN = log(2), summary = TRUE, fileBase = "result")}
\arguments{
  \item{phy}{a phylogenetic tree of class 'phylo'}
  \item{dat}{a named vector of continuous trait values, associated with each species in \code{phy}; see \code{\link[geiger]{name.check}}}
  \item{SE}{a named vector of standard errors for each trait value; applied to all trait values if given a single value}
  \item{ngen}{number of sampling generations}
  \item{sample.freq}{frequency with which Markov samples are retained (e.g., \code{sample.freq=10} retains every tenth sample in the chain)}
  \item{prob.lambda}{proportion of proposals that split or merge rates, thereby altering complexity of the evolutionary model}
  \item{prob.sigsqB}{proportion of proposals in which the root state is updated}
  \item{prob.sigsqJ}{a constraint of model complexity (e.g., \code{constrainK=2} restricts sampling to models with two independent rates)}
  \item{prop.width}{if \code{NULL}, the proposal width is calibrated; otherwise, proposal width is constrained to the supplied value}
  \item{lambdaN}{the shape parameter for the Poisson prior-distribution on number of distinct rates}
  \item{summary}{a logical switch that determines whether any output is printed to the console after the analysis has finished}
  \item{fileBase}{a string used to uniquely identify the base directory to which Markov samples are written}
}

\details{
In the two primary objects that store posterior estimates of \code{shifts} and \code{rates} (see \bold{Value}), branches are labeled in accordance with the numeric identifiers found in \code{phy$edge}, which represents a table of ancestor-descendant
 relationship (see \code{\link[ape]{read.tree}}.  In the \bold{Examples} below do \bold{not} expect convergence with such short chains!
}
\value{
If \code{summary=TRUE} the global rate of acceptance and calibrated proposal width is returned as a list.  After a run has completed, the rates of acceptance for various proposal mechanisms are printed, along with 
the (calibrated) proposal width used for chain sampling. 

Posterior results are written to several files within a base directory, the contents of which are as follows:
  \item{logfile}{a logfile including the following for each Markov sample; the generations at which samples were retained; the generalized model of trait evolution used for inference (e.g., Brownian motion or 'BM"); 
  the mean treewide rate in the sampled generation; the number of independent rates; the root state; and the likelihood of the model at each sample.  Do not be alarmed if it seems to take some time before values 
  begin to appear in the \code{logfile}: samples are not saved until after the calibration period has terminated. Inspection of chain mixing can be facilitated by the R-package \pkg{coda} or
   by the Java application, Tracer (\url{http://tree.bio.ed.ac.uk/software/tracer/}). }
  \item{errors}{a logfile to which errors are recorded, if any are generated in the Markov chain}
  \item{posteriorsamples}{a compressed .rda file, storing estimates from the Markov chain for relevant parameters; contents of this file can be loaded 
  with \code{load("posteriorsamples.rda")}. This list object, whose referent is also \code{posteriorsamples} when in R-memory, includes two dataframes: 
  \code{rates} comprises the relative-rate estimates for each branch in the phylogeny; \code{rate.shifts} stores the branches that were chosen for 
  rate shifts in each retained sample of the Markov chain (\code{1} signifies a shift; 
  \code{0} implies a branch that was not sampled for a rate shift). Column names of both components of \code{posteriorsamples} are from the numeric vector of 
  branch labels (see \bold{Details}). For further information on .Rd and .rda files, see \code{\link[base]{save}} and \code{\link[base]{load}} as well as the 
  \bold{Examples} below. }
}

\author{JM Eastman}
\examples{
#############
## generate data
tmp=phenogram(rcoal(200), "levy"); phy=tmp$phy; dat=tmp$dat
mcmc.levy(phy,dat,prop.width=0.5)
mcmc.levy(phy,dat)
}
