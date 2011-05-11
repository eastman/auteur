require(auteur)
bytens=function(x, num) sapply(0:num, function(y)x*(10^y))
vcv2 <- function(tree){
    res = cophenetic(tree) / 2
    max(res) - res
}

n=bytens(20,2)
reps=5
out=array(dim=c(length(n)*reps,5))
for(i in n){
	for(j in 1:reps){
		cat(paste("working on rep", j, "of size", i, "\n",sep=" "))
		phy=rcoal(i)
		ii=which(n==i)
		sA=system.time(vcv(phy)->vA)[3] 
		sR=system.time(vmat(phy)->vR)[3] 
		s2=system.time(vcv2(phy)->v2)[3]
		match=all(vA==vR&vA==v2) 
		out[(ii-1)*reps+j,]=c(match, sA,sR,s2,Ntip(phy))
		}
		
		}
out=data.frame(out)
out=cbind(out, out[,2]/out[,3],out[,2]/out[,4])
names(out)=c("match","ape","rjmcmc","klaus","treesize","speedupRJMCMC","speedupKLAUS")
save(out,file="checkVMAT.rda")
