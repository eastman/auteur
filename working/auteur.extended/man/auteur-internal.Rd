\name{auteur-internal}
\alias{add.transparency}
\alias{adjustrate}
\alias{adjustvalue}
\alias{assess.lnR}
\alias{bm.lik.fn}
\alias{branchcol.plot} 
\alias{choose.one}
\alias{cleanup.files}
\alias{collect.nodes}
\alias{determine.accepted.proposal}
\alias{edgelabels.rj}
\alias{find.relatives}
\alias{fit.continuous}
\alias{generate.error.message}
\alias{generate.log}
\alias{generate.starting.point}
\alias{get.ancestor.of.node}
\alias{get.ancestors.of.node}
\alias{get.desc.of.node}
\alias{get.descendants.of.node}
\alias{hdr}
\alias{infer.y}
\alias{is.root}
\alias{nearest.pair}
\alias{parlogger}
\alias{pathlengths}
\alias{prepare.data}
\alias{prunelastsplit}
\alias{ptpois}
\alias{rtpois}
\alias{splitormerge}
\alias{splitrate}
\alias{splitvalue}
\alias{summarize.run}
\alias{table.print}
\alias{tune.rate}
\alias{tune.value}
\alias{updatevcv}
\alias{ultrametricize}
\alias{withinrange}

\title{internal auteur functions}
\description{\code{auteur-internal} functions are either not typically called by the user or are currently undocumented}

\details{
This is an internal \pkg{auteur} function, not intended to be called directly by the user.
}

\note{The likelihood expression used in \code{bm.lik.fn} owes to Felsenstein (1973); see also O'meara et al. (2006)}

\references{
FELSENSTEIN J. 1973. Maximum likelihood estimation of evolutionary trees from continuous characters. American Journal of Human Genetics 25:471-492.
O'MEARA BC, C ANE, MJ SANDERSON, and PC WAINWRIGHT. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60:922-933.
}
