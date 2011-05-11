#rjmcmc run diagnosis, generating plots to compare posterior and prior densities of K-rate models and conduct Bayes factor comparisons between models differing in complexity
#author: JM EASTMAN 2010

tracer <-
function(base.log, lambdaK=log(2), Bayes.factor=NULL, burnin=0.25, col.line=1, pdf=TRUE, factor="M",...) {
	d=read.table(base.log,header=TRUE)
	tag.tmp=strsplit(base.log,".",fixed=TRUE)[[1]]
	tag=paste(tag.tmp[(tag.tmp != "rjTrait" & tag.tmp != "log")],collapse=".")
	
	if(pdf) pdf(paste(tag,"pdf",sep="."))
		
	# lnL trace
	if(any(names(d)=="lnL")){
		if(!is.null(factor)) {
			if(factor=="k" | factor=="K") {
				plot(d$gen/1000, d$lnL, main="likelihood trace", ylab="log-likelihood", xlab="generation (k)", bty="n", cex=0.5)
			} else {
				plot(d$gen/1000000, d$lnL, main="likelihood trace", ylab="log-likelihood", xlab="generation (M)", bty="n", cex=0.5)
			}
		} else {
			plot(d$gen, d$lnL, main="likelihood trace", ylab="log-likelihood", xlab="generation", bty="n", cex=0.5)	
		}
	}
	
	# mean rate density
	if(any(names(d)=="mean.rate")) dmr=as.data.frame(d$mean.rate) 
	if(length(col.line)!=ncol(dmr)) {
		col.line=gray.colors(ncol(dmr)+1)
		col.line=col.line[-length(col.line)]
	}
	if(!pdf) dev.new()
	if(all(dmr>0)) {
		trace.plot(dmr, col=col.line, xlab="rate", hpd=NULL, legend.control=list(plot=FALSE), main="posterior sample of mean rates", log="x", ...)	
	} else {
		trace.plot(dmr, col=col.line, xlab="rate", hpd=NULL, legend.control=list(plot=FALSE), main="posterior sample of mean rates", ...)	
	}
	
	# K-rate model prior and posterior probabilities; Bayes factor comparison
	if(any(names(d)=="rates")) {
		bb=max(1, round(burnin*nrow(d)))
		d=d[bb:nrow(d),]
		prior=rpois(nrow(d), lambdaK)+1
		tt.prior=table(prior)
		tt.post=table(d$rates)
		df=array(dim=c(2,max(c(prior,d$rates))))
		for(pr in 1:length(tt.prior)) {
			df[1,(as.numeric(names(tt.prior)[pr]))]=tt.prior[pr]
		}
		for(pp in 1:length(tt.post)) {
			df[2,(as.numeric(names(tt.post)[pp]))]=tt.post[pp]
		}
		df[is.na(df)]=0
		dimnames(df)[[2]]=seq(1,ncol(df))
		if(!pdf) dev.new()
		df=apply(df,2,function(x)x/sum(df/2))
		
		# posterior-prior plot of K-rate models
		barplot(df,col=c("whitesmoke","black"),xlab="rates",ylab="probability",ylim=c(0,1),beside=TRUE,main="probabilities of k-rate models")
		legend(x="topright", c("prior","posterior"), fill=c("whitesmoke","black"), bty="n", cex=0.75)
		
		# Bayes factor comparison
		if(!is.null(Bayes.factor)) {
			is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
			resample <- function(a=c(0,0), b=vector(), resamples=1000, frequency=0.5, multi=FALSE) {
				if(!multi) {
					xx=lapply(1:resamples, function(x) {
							  rr=sample(b, size=round(frequency*length(b)), replace=TRUE)
							  return(sum(rr==a[2])/sum(rr==a[1]))
							  }
							  )
				} else {
					xx=lapply(1:resamples, function(x) {
							  rr=sample(b, size=round(frequency*length(b)), replace=TRUE)
							  return(sum(rr>1)/sum(rr==1))
							  }
							  )
				}
				return(unlist(xx))
			}
			if(length(Bayes.factor)==2 ) {
				if(!pdf) dev.new()
				if(Bayes.factor[2]!="multi") {
					BF.tmp=sapply(Bayes.factor, function(x)if(sum(d$rates==x)>0)return(TRUE)else return(FALSE))
					if(all(BF.tmp)) {
						BF=sum(d$rates==Bayes.factor[2])/sum(d$rates==Bayes.factor[1])
						resampled.BF=resample(Bayes.factor, d$rates)
						trace.plot(resampled.BF[is.finite(resampled.BF)], legend.control=list(plot=FALSE), main=paste("Bayes factor: ", paste(paste(rev(Bayes.factor), "-", sep=""), collapse=" versus "), "rate model of evolution", sep=""), xlab="Bayes factor")
						lines(data.frame(c(BF,BF),c(0,max(density(resampled.BF)$y)), lty=4))
					} else {
						BF=paste("Model(s)", paste(Bayes.factor[which(!BF.tmp)], collapse=" and "), "appear(s) to be missing from the posterior sample.", sep=" ")
						resampled.BF=NA
					}
				} else {
					BF.tmp=c((sum(d$rates==1)>0), (sum(d$rates>1)>0))
					if(all(BF.tmp)) {
						BF=sum(d$rates>1)/sum(d$rates==1)
						resampled.BF=resample(Bayes.factor, d$rates, multi=TRUE)
						trace.plot(resampled.BF[is.finite(resampled.BF)], legend.control=list(plot=FALSE), main=paste("Bayes factor: ", paste(paste(c("multi","global"), "-", sep=""), collapse=" versus "), "rate model of evolution", sep=""), xlab="Bayes factor")
						lines(data.frame(c(BF,BF),c(0,max(density(resampled.BF)$y)), lty=4))
					} else {
						BF=paste("Model(s)", paste(Bayes.factor[which(!BF.tmp)], collapse=" and "), "appear(s) to be missing from the posterior sample.", sep=" ")
						resampled.BF=NA
					}
				}
				
			}
		}
	}

	if(pdf) dev.off()	
	if(!is.null(Bayes.factor)) return(list(BF=BF, resampled.BF=resampled.BF))
}

