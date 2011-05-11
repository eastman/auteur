#general phylogenetic plotting utility, which is a modification of ape:::edgelabels, plotting edge symbols at the temporal beginning of the branch
#author: E PARADIS 2009 and JM EASTMAN 2010
#note: may not be trustworthy where lastPP$type is not "phylogram"

edgelabels.rj <-
function (text, edge, adj = c(0.5, 0.5), frame = "rect", pch = NULL, 
thermo = NULL, pie = NULL, piecol = NULL, col = "black", 
bg = "lightgreen", horiz = FALSE, width = NULL, height = NULL, 
...) 
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(edge)) {
        sel <- 1:dim(lastPP$edge)[1]
        subedge <- lastPP$edge
    }
    else {
        sel <- edge
        subedge <- lastPP$edge[sel, , drop = FALSE]
    }
    if (lastPP$type == "phylogram") {
        if (lastPP$direction %in% c("rightwards", "leftwards")) {
            XX <- (lastPP$xx[subedge[, 1]])
            YY <- lastPP$yy[subedge[, 2]]
        }
        else {
            XX <- lastPP$xx[subedge[, 2]]
            YY <- (lastPP$yy[subedge[, 1]])
        }
    }
    else {
        XX <- (lastPP$xx[subedge[, 2]])
        YY <- (lastPP$yy[subedge[, 2]])
    }
	if(missing(text)) text=lastPP$edge[,2]

    ape:::BOTHlabels(text, sel, XX, YY, adj, frame, pch, thermo, pie, 
					 piecol, col, bg, horiz, width, height, ...)
}

