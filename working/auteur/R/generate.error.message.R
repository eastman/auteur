#logging utility used by rjmcmc
#author: JM EASTMAN 2010

generate.error.message <-
function(ntip, i, mod.cur, mod.new, lnL, lnPrior, lnHastings, curCats, newCats, errorLog) {
	if(!file.exists(file=errorLog)) write(paste("nodes", "gen", "cur.lnL", "new.lnL", "lnL", "lnPrior", "lnHastings", "cur.catg", "new.catg", sep="\t"), file=errorLog)

	write(paste(2*ntip-2, i, sprintf("%.3f", mod.cur$lnL), sprintf("%.3f", mod.new$lnL), sprintf("%.3f", lnL), sprintf("%.3f", lnPrior), sprintf("%.3f", lnHastings), sum(curCats), sum(newCats), sep="\t"),  file=errorLog, append=TRUE)
	
}

