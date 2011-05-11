#rjmcmc run diagnosis, generating tallies of proposed and accepted updates by class of proposal mechanism
#author: JM EASTMAN 2010

summarize.run <-
function(cOK, endparms, cur.model) {
	if(cur.model=="BM") {
		prop.names<-c("rate.adj", "rate.tune", "rate.catg", "root")
	} else {
		prop.names<-c("rate", "root", "alpha", "regimes", "theta")
	}
	df=data.frame(cbind(proposed=endparms, accepted=cOK, adoptrate=cOK/endparms))
	rownames(df)=prop.names
	cat("\n\n",rep(" ",10),toupper(" sampling summary"),"\n")

	if(any(is.na(df))) df[is.na(df)]=0
	table.print(df, digits=c(0,0,4), buffer=6)
	
	cat("\n\n")

}

