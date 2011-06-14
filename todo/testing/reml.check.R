require(auteur)

orig.bm.reml.fn <- 
function(phy, dat, sigsq) { # from LJ HARMON and LJ REVELL
	ic<-pic(dat, phy)
	reml<-dnorm(ic, mean=0, sd=sqrt(sigsq), log=TRUE)
	return(sum(reml))
}


#new.bm.reml.fn <- 
#function(phy, dat, rates, orig.var, ic) {
#	new.tre=updatetree(phy, rates)
#	new.var=pic(dat, new.tre, scale=FALSE, var=TRUE)[,"variance"]
#	reml=dnorm(ic, mean=0, sd=sqrt(new.var/orig.var), log=TRUE)
#	return(sum(reml))
#}

updatetree <-
function(ape.tre, new.rates){
	ape.tre$edge.length=ape.tre$edge.length*new.rates
	return(ape.tre)
}


# simulate data
phy=rcoal(20)
fast.phy=phy
fast.phy$edge.length=fast.phy$edge.length*10
dat=sim.char(fast.phy,model="brownian",model.matrix=as.matrix(0.5))[,,1]

rr=seq(.1, 100, by=0.5)

# standard ML
yl=sapply(rr, function(x) {v=updatevcv(phy, rep(x, length(phy$edge.length))); bm.lik.fn(root=0, dat=dat, vcv=v, SE=NULL)$lnL})
plot(rr,yl, ylim=c(max(yl-20), max(yl)), main=paste("ML relrate: ", rr[yl==max(yl)], sep=""))

# single rate REML -- works fine [ specify rate via sigmasq ]
yr=sapply(rr, function(x) {orig.bm.reml.fn(phy, dat, sigsq=x)})
dev.new()
plot(rr,yr, ylim=c(max(yr-20), max(yr)), main=paste("REML relrate: ", rr[yr==max(yr)], sep=""))

unscl=pic(dat,phy,scaled=FALSE)
orig.var=pic_variance(phy)

# generalized REML -- working [ tree transformation: find ratio of transformed expected variances to original expected variances in order to find new rate]
yg=sapply(rr, function(x) {tt=phy; tt$edge.length=tt$edge.length*x; new.var=pic(dat, tt, scale=FALSE, var=TRUE)[,2]; sum(dnorm(unscl, mean=0, sd=sqrt(new.var), log=TRUE))})
dev.new()
plot(rr,yg, ylim=c(max(yg-20), max(yg)), main=paste("gREML relrate: ", rr[yg==max(yg)], sep=""))

# generalized REML -- working [ tree transformation: find ratio of transformed expected variances to original expected variances in order to find new rate]
yt=sapply(rr, function(x) bm.reml.fn(reorder(phy,"pruningwise"), rep(x, length(phy$edge.length)), unscl)$lnL)
dev.new()
plot(rr,yt, ylim=c(max(yt-20), max(yt)), main=paste("tREML relrate: ", rr[yt==max(yt)], sep=""))

