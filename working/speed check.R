#####
packs=c("rjmcmc.traits_0.11.0106","auteur_0.11.0119")
reps=5
n=c(10,100,1000)
out=array(dim=c(length(packs)*reps*length(n),2))
library(rjmcmc.traits)
taxon="urodela"
ngen=2000

for(i in packs) {
	if(i=="auteur_0.11.0119") detach("package:rjmcmc.traits", unload=TRUE)
	install.packages(paste(paste(i,"tar.gz",sep="."),sep="/"),repos=NULL,type="source")	
	if(i=="auteur_0.11.0119") library(rjmcmc.traits) else library(auteur)	
	for(nn in n) {
		for(j in 1:reps) {
			ind=which(packs==i)
			phy=rcoal(nn)
			dat=rTraitCont(phy)
			out[(reps*(ind-1)+j)*(which(nn==n)),]=c(system.time(rjmcmc.bm(phy, dat, ngen=ngen, sample.freq=50, prob.mergesplit = 0.10, prob.root=0.005, jumpsize=1, internal.only=FALSE, fileBase = paste(taxon,i,sep=".")))[3],i)
		}
		
		
	}
	
}

out=data.frame(out)
names(out)=c("time","package")
print(out)
save(out,file="speed.rda")
