require(auteur)
byN=function(x, num, N) sapply(0:num, function(y)x*(N^y))
n=byN(20, 8, 2)
reps=3
out=array(dim=c(length(n)*reps,2))
for(i in n){
	for(j in 1:reps){
		cat(paste("working on rep", j, "of size", i, "\n",sep=" "))
		phy=rcoal(i)
		ii=which(n==i)
		out[(ii-1)*reps+j,]=c(i, system.time(vmat(phy))[3])
	}		
}
out=data.frame(out)
names(out)=c("treesize","time")
save(out,file="timeVMAT.rda")
