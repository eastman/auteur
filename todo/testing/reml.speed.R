require(auteur)
reps=2^(6:12)
dat=data.frame(N=reps, reml=NA, orig=NA)
for(i in 1:nrow(dat)){
	cat(paste("WORKING ON SET: ", dat$N[i], "\n", sep=""))
	phy=rcoal(dat$N[i])
	trt=rTraitCont(phy)
	ptre=reorder(phy, "pruningwise")
	ic=pic(trt, ptre, scaled=FALSE) 
	dat$reml[i]=system.time(bm.reml.fn(ptre, rep(1, length(phy$edge.length)), ic))[3]
	dat$orig[i]=system.time(bm.lik.fn(0, trt, vmat(phy), 0))[3]		
}

dat$reml[dat$reml==0]=0.0001
dat$orig/dat$reml->rel
plot(dat$N,rel/1000, log="x", xlab="treesize", ylab="speedup (x 1000)")
