inverseA.HadfieldNakagawa<-function(phy=NULL){
	require(Matrix)
	
## CONVERT TO C++
# currently vmat is faster for vcv computation (but deals with NxN matrix, where N is tips)
# using solve on the vmat calculated vcv matrix is slower than inverseA of MCMCglmm
# need to combine the efficiency of building the matrix in C++ (here vectorized) and the efficiency of that introduced by Hadfield and Nakagawa (2010), where they deal with
# the inverse of S (nxn) rather than the inverse of A (NxN), where n is all nodes
# further efficiency seems to be lent by dealing with sparse matrices explicitly in this function
	
    roots<-setdiff(phy$edge[,1], phy$edge[,2])  
	
# reordering done to allow vectorized computation (I think)
    reorder<-order((phy$edge[,1]%in%roots==FALSE)*phy$edge[,2]+10000000*(phy$edge[,2]<=length(phy$tip.label)))
	
    id<-phy$edge[,2][reorder]
    tips<-which(id<=length(phy$tip.label))
	
    node.names<-id
    node.names[tips]<-phy$tip.label
	
    if(is.null(phy$node.label)){
		phy<-makeNodeLabel(phy) 
    }
	
    node.names[-tips]<-phy$node.label[-1]
	
    edges<-phy$edge.length[reorder]
    if(any(edges<1e-16)){stop("some phylogeny edge.lengths are zero")}
	
    dam<-match(phy$edge[,1][reorder], id)
    id<-match(id, id)
	
#   phy<-cbind(node.names, node.names[dam], rep(NA, length(node.names)))
	
    dam[which(is.na(dam))]<--998
	
	
    Ainv<-Matrix(0, length(dam), length(dam))
    off<-tapply(id, dam,function(x){x})
    off<-off[-which(names(off)=="-998")]
    off<-lapply(off, function(x){sum(1/edges[x], na.rm=T)})
    
## TO C++
	diag(Ainv)<-as.numeric(1/edges)
    diag(Ainv)[as.numeric(names(off))]<-diag(Ainv)[as.numeric(names(off))]+unlist(off)
    Ainv[(dam[which(dam!="-998")]-1)*dim(Ainv)[1]+id[which(dam!="-998")]]<-(-1/edges[id[which(dam!="-998")]])
    Ainv[(id[which(dam!="-998")]-1)*dim(Ainv)[1]+dam[which(dam!="-998")]]<-(-1/edges[id[which(dam!="-998")]])
## END to C++
	
    Ainv.tips<-as(Ainv[,tips][tips,]-Ainv[,-tips][tips,]%*%as(solve(as.matrix(Ainv[,-tips][-tips,])), "sparseMatrix")%*%Ainv[,tips][-tips,],"matrix")	
	mm=match(node.names[tips],phy$tip.label)
	Ainv.tips=Ainv.tips[mm,mm]
	
	
 
   
  return(Ainv.tips)
}

