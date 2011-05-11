function (phy, data, ngen = 10000, probMergeSplit = 0.01, rootp = 0.1, 
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
    }
    else {
        f1 <- ""
        f2 <- ""
    }
    cat("", file = f1)
    cat("", file = f2)
    for (i in 1:ngen) {
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
	## try added in next line
        lnlNew <- try(fitModel(td$phy, td$data[, 1], newRate, newRootMean), silent = F)
        ## if the proposal generates an error, treats the likelihood as zero (i.e., just keeps the original) and logs the error
        ## this should be changed so that if the proposal generates an error, a new proposal is offered in the same generation
        if (class(lnlNew) == "try-error" || is.na(lnlNew)) {
            errorNum <- errorNum + 1
            cat(c('Error number', errorNum, 'occurred in generation', i), file = errorLog, append = T, sep = ' ')
            if(is.na(lnlNew)) {cat(c("Likelihood calc NaN error", '\n'), file = errorLog, append = T, sep = ' ')}
	      else cat(c(lnlNew, '\n'), file = errorLog, append = T, sep = ' ')
            cat(c('Proposed rate and mean:', '\n'), file = errorLog, append = T, sep = ' ')
            cat(c(newRate, '\n'), file = errorLog, append = T, sep = ' ')
            cat(c(newRootMean, '\n\n'), file = errorLog, append = T, sep = ' ')
            cat(c('Previous rate and mean:', '\n'), file = errorLog, append = T, sep = ' ')
            cat(c(currRate, '\n'), file = errorLog, append = T, sep = ' ')
            cat(c(currRootMean, '\n'), file = errorLog, append = T, sep = ' ')
            cat(c('--------------------------------------------', '\n\n'), file = errorLog, append = T, sep = ' ')
            cat(c(i, '\t', length(unique(newRateCat)),'\t', min(newRate), '\t', max(newRate), '\t', round(errorNum / i, 2), '\n'), file = paste(errorLog, 'stats.txt', sep = ''), append = T, sep = '')
            res$lnlprof[i] <- lnlCurr
        }
        else {
            lnLikelihoodRatio = lnlNew - lnlCurr
            message(paste('lnlNew - lnlCurr =', lnlNew, '-', lnlCurr, '=', lnLikelihoodRatio))
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
            }
            else {
                res$lnlprof[i] <- lnlCurr
            }
        }
        if (i%%sampleFreq == 0) {
            res$rate <- cbind(res$rate, currRate)
            res$rateCat <- cbind(res$rateCat, currRateCat)
            cat(file = f1, lnlCurr, currRateCat, "\n", append = T, 
                sep = "\t")
            cat(file = f2, currRate, "\n", append = T, sep = "\t")
        }
        if (i%%100 == 0) {
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
