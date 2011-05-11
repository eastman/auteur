#general phylogenetic utility for selecting one element of a list of values, each of which can be associated with the edges of a phylogeny (phy$edge[,2])
#author: JM EASTMAN 2010

choose.one <-
function(cur.delta, phy=NULL, internal.only=FALSE)
# updated 12.16.2010
{
	bb=cur.delta
	if(internal.only) {
		while(1) {
			s=sample(1:length(bb),1)
#			print(bb[s])
			names(bb)=phy$edge[,2]
			if(as.numeric(names(bb[s]))>Ntip(phy)) break()
		}
	} else {
		s=sample(1:length(bb),1)
	}
	bb[s]=1-bb[s]
	bb
}

