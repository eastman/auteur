brownMcmc<-function (phy, data, ngen = 10000, probMergeSplit = 0.01, rootp = 0.1, 
    sampleFreq = 1000, ylim = c(-200, 0), fileBase = "") 
## Version 5: only difference is in lines 157 -- 165, using 'try' to catch bad vcv matrices (AH and LH, 22 apr 2008)
{
    ## Added 18 june 08, AH
      errorNum <- 0 
      errorLog <- paste(fileBase, 'errors.log', sep = '')
      cat(c('Error log for run beginning', date(), '\n'), file = errorLog, append = F, sep = ' ')
      cat(c('--------------------------------------------------------', '\n\n'), file = errorLog, append = T)
      cat(c('generation', '\t', 'numCat', '\t', 'rateMin', '\t', 'rateMax', '\t', 'cumulativePercentage', '\n'), file = paste(errorLog, 'stats.txt', sep = ''), append = F, sep = '')
    td <- treedata(phy, data, sort = T)
    nn <- length(td$phy$edge.length)
    currRateCat <- numeric(nn)
    currRateCat[] <- 1
    currRate <- numeric(nn)
    currRate[] <- 0.001
    currRootMean <- mean(td$data[, 1])
    l <- numeric(ngen)
    l[] <- NA
    nMergeProp = 0
    nSplitProp = 0
    nRateProp = 0
    nMergeOK = 0
    nSplitOK = 0
    nRateOK = 0
    res <- list(rate = currRate, rateCat = currRateCat, lnlprof = l)
    if (fileBase != "") {
        f1 <- paste(fileBase, "RateCat.txt", sep = "")
        f2 <- paste(fileBase, "Rates.txt", sep = "")
    } else {
        f1 <- ""
        f2 <- ""
    }
    cat("", file = f1)
    cat("", file = f2)
    for (i in 1:ngen) {
		print(max(currRateCat))
        lnLikelihoodRatio = 0
        lnProposalRatio = 0
        lnPriorRatio = 0
        merge = FALSE
        split = FALSE
        if (runif(1) < (2 * probMergeSplit)) {
            if (max(currRateCat) == 1) {
                split = TRUE
            }
            else if (max(currRateCat) == nn) {
                merge = TRUE
            }
            else if (runif(1) < 0.5) {
                split = TRUE
            }
            else merge = TRUE
        }
        if (merge) 
            nMergeProp = nMergeProp + 1
        if (split) 
            nSplitProp = nSplitProp + 1
        if (!merge & !split) 
            nRateProp = nRateProp + 1
        newRateCat <- currRateCat
        newRate <- currRate
        newRootMean <- currRootMean
        if (!merge & !split) {
            if (runif(1) < rootp) {
                rch <- runif(1, min = -0.5, max = 0.5)
                newRootMean <- newRootMean + rch
            }
            else {
                ncat <- max(currRateCat)
                rcat <- round(runif(1) * (ncat) + 0.5)
                r <- log(currRate[currRateCat == rcat][1])
                rch <- runif(1, min = -2, max = 2)
                nr <- r + rch
                newRate[currRateCat == rcat] <- exp(nr)
            }
        }
        if (merge) {
            ncat <- max(currRateCat)
            mergeCats <- sample(1:ncat)[1:2]
            newRateCat[newRateCat == mergeCats[1]] <- mergeCats[2]
            newRateCat <- fixRateCat(newRateCat)
            or1 <- currRate[currRateCat == mergeCats[1]][1]
            or2 <- currRate[currRateCat == mergeCats[2]][1]
            n1 <- sum(currRateCat == mergeCats[1])
            n2 <- sum(currRateCat == mergeCats[2])
            nr <- (n1 * or1 + n2 * or2)/(n1 + n2)
            newRate[currRateCat == mergeCats[1] | currRateCat == 
                mergeCats[2]] <- nr
            numSplittableBeforeMerge <- sum(table(currRateCat) > 
                1)
            numSplittableAfterMerge <- sum(table(newRateCat) > 
                1)
            if (max(currRateCat) == nn) {
                probMerge = 1
            }
            else {
                probMerge = 0.5
            }
            if (max(newRateCat) == 1) {
                probSplit = 1
            }
            else {
                probSplit = 0.5
            }
            factor = probSplit/probMerge
            lnProposalRatio = log(factor) + log(nn) + log(ncat) - 
                log(numSplittableAfterMerge) - log(2^(n1 + n2) - 
                2) - log(nr/mean(currRate)) - log(n1 + n2)
            lnPriorRatio = log(Stirling2(nn, ncat)) - log(Stirling2(nn, 
                ncat - 1))
        }
        if (split) {
            ncat <- max(currRateCat)
            while (1) {
                rcat <- round(runif(1) * ncat + 0.5)
                nc <- sum(currRateCat == rcat)
                if (nc > 1) 
                  break
            }
            while (1) {
                new <- round(runif(nc) * 2 + 0.5)
                if (length(table(new)) == 2) 
                  break
            }
            newRateCat[currRateCat == rcat][new == 2] <- ncat + 
                1
            or <- currRate[currRateCat == rcat][1]
            n1 <- sum(new == 1)
            n2 <- sum(new == 2)
            u <- runif(1, min = -0.5 * n1 * or, max = 0.5 * n2 * 
                or)
            nr1 <- or + u/n1
            nr2 <- or - u/n2
            newRate[newRateCat == rcat] <- nr1
            newRate[newRateCat == ncat + 1] <- nr2
            newRateCat <- fixRateCat(newRateCat)
            numSplittableBeforeSplit <- sum(table(currRateCat) > 
                1)
            numSplittableAfterSplit <- sum(table(newRateCat) > 
                1)
            if (max(currRateCat) == 1) {
                probSplit = 1
            }
            else {
                probSplit = 0.5
            }
            if (max(newRateCat) == nn) {
                probMerge = 1
            }
            else {
                probMerge = 0.5
            }
            factor = probMerge/probSplit
            lnProposalRatio = log(factor) + log(numSplittableBeforeSplit) + 
                log(2^(n1 + n2) - 2) + log(or/mean(currRate)) + 
                log(n1 + n2) - log(nn) - log(ncat + 1)
            lnPriorRatio = log(Stirling2(nn, ncat)) - log(Stirling2(nn, 
                ncat + 1))
        }
        lnlCurr <- fitModel(td$phy, td$data[, 1], currRate, currRootMean)
        lnlNew <-  fitModel(td$phy, td$data[, 1], newRate, newRootMean)

            lnLikelihoodRatio = lnlNew - lnlCurr
            heat = 1
            lnR = heat * lnLikelihoodRatio + heat * lnPriorRatio + 
                lnProposalRatio
    
            if (lnR > 0) {
                r = 1
            }
            else if (lnR < -100) {
                r = 0
            }
            else {
                r = exp(lnR)
            }
            #message(paste('lnlNew - lnlCurr =', lnlNew, '-', lnlCurr, '=', lnLikelihoodRatio, "lnR", lnR, "r", r))
            
            if (runif(1) <= r) {
                currRate <- newRate
                currRateCat <- newRateCat
                currRootMean <- newRootMean
                res$lnlprof[i] <- lnlNew
                if (merge) 
                  nMergeOK = nMergeOK + 1
                if (split) 
                  nSplitOK = nSplitOK + 1
                if (!merge & !split) 
                  nRateOK = nRateOK + 1
                #message("YES")  
            }
            else {
                res$lnlprof[i] <- lnlCurr
                #message("NO")  

            }
        
        if (i%%sampleFreq == 0) {
            res$rate <- cbind(res$rate, currRate)
            res$rateCat <- cbind(res$rateCat, currRateCat)
            cat(file = f1, lnlCurr, currRateCat, "\n", append = T, 
                sep = "\t")
            cat(file = f2, currRate, "\n", append = T, sep = "\t")
        }
        if (i%%1000 == 0) {
            cat("Gen ", i, "completed\n")
        }
    }
    res$nMergeOK = nMergeOK
    res$nMergeProp = nMergeProp
    res$nSplitOK = nSplitOK
    res$nSplitProp = nSplitProp
    res$nRateOK = nRateOK
    res$nRateProp = nRateProp
    res
}


fixRateCat<-function(x)
{
	f<-factor(x)
	n<-nlevels(f)
	o<-numeric(n)
	for(i in 1:n) o[i]<-which(f==levels(f)[i])[1]
	nf<-ordered(f, levels(f)[order(o)])	
	nx<-as.numeric(nf)
	nx
}

fitModel<-function(phy, data, model, rootMean)
{
	phy$edge.length=phy$edge.length*model
	vcv<-vcv.phylo(phy)
	mm<-rep(rootMean, nrow(vcv))
	return(lnldmv(data, mm, vcv))
	
	}
	
Stirling2 <- function(n,m)
{
    ## Purpose:  Stirling Numbers of the 2-nd kind
    ## 		S^{(m)}_n = number of ways of partitioning a set of
    ##                      $n$ elements into $m$ non-empty subsets
    ## Author: Martin Maechler, Date:  May 28 1992, 23:42
    ## ----------------------------------------------------------------
    ## Abramowitz/Stegun: 24,1,4 (p. 824-5 ; Table 24.4, p.835)
    ## Closed Form : p.824 "C."
    ## ----------------------------------------------------------------

    if (0 > m || m > n) stop("'m' must be in 0..n !")
    k <- 0:m
    sig <- rep(c(1,-1)*(-1)^m, length= m+1)# 1 for m=0; -1 1 (m=1)
    ## The following gives rounding errors for (25,5) :
    ## r <- sum( sig * k^n /(gamma(k+1)*gamma(m+1-k)) )
    ga <- gamma(k+1)
    round(sum( sig * k^n /(ga * rev(ga))))
}

lnldmv<-function(y, mean, sigma) {
	num<- -t(y-mean)%*%solve(sigma)%*%(y-mean)/2
	den<-sqrt((2*pi)^length(y)*det(sigma))
	return((num-log(den))[1,1])
	}	
