\name{auteur-internal}
\alias{add.transparency}
\alias{adjustjump}
\alias{adjustrate}
\alias{adjustvalue}
\alias{assess.lnR}
\alias{assigndescendants}
\alias{bm.lik.fn}
\alias{branchcol.plot} 
\alias{choose.one}
\alias{cleanup.files}
\alias{collect.nodes}
\alias{determine.accepted.proposal}
\alias{edgelabels.rj}
\alias{find.relatives}
\alias{fit.continuous}
\alias{fit.hansen}
\alias{generate.error.message}
\alias{generate.log}
\alias{generate.starting.point}
\alias{get.ancestor.of.node}
\alias{get.ancestors}
\alias{get.ancestors.of.node}
\alias{get.desc.of.node}
\alias{get.descendants.of.node}
\alias{get.epochs}
\alias{get.times}
\alias{hdr}
\alias{infer.y}
\alias{is.root}
\alias{jumps.edgewise}
\alias{levyevolver}
\alias{nearest.pair}
\alias{parlogger}
\alias{pathlengths}
\alias{prepare.data}
\alias{priorratio.exp}
\alias{prunelastsplit}
\alias{ptpois}
\alias{ou.lik.fn}
\alias{oumat}
\alias{ouweights}
\alias{reshape.localvalues}
\alias{rtpois}
\alias{splitormerge}
\alias{splitrate}
\alias{splitvalue}
\alias{summarize.run}
\alias{table.print}
\alias{tree.slide}
\alias{tree.slide.LEVY}
\alias{tune.rate}
\alias{tune.value}
\alias{updatevcv}
\alias{updateJvcv}
\alias{updateLvcv}
\alias{ultrametricize}
\alias{withinrange}

\title{internal auteur functions}
\description{\code{auteur-internal} functions are either not typically called by the user or are currently undocumented}

\details{
This is an internal \pkg{auteur} function, not intended to be called directly by the user.
}

\note{The likelihood expression used in \code{bm.lik.fn} owes to Felsenstein (1973), and see also O'meara et al. (2006). That used within \code{ou.lik.fn} owes to Butler and King (2004).}

\references{
\itemize{
	\item{BUTLER MA and AA KING. 2004. Phylogenetic comparative analysis: A modeling approach for adaptive evolution. American Naturalist 164:683-695.}
	\item{FELSENSTEIN J. 1973. Maximum likelihood estimation of evolutionary trees from continuous characters. American Journal of Human Genetics 25:471-492.}
	\item{O'MEARA BC, C ANE, MJ SANDERSON, and PC WAINWRIGHT. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60:922-933.}
}
}
