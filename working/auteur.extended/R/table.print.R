#general printing utility for ensuring equal numbers of characters within columns and defining spacing between columns
#author: JM EASTMAN 2010
#note: works only for numeric dataframes

table.print=function(df,digits=4,buffer=5){
	if(length(buffer) != ncol(df) | length(buffer)==1) buffer=rep(buffer[1],ncol(df))
	if(length(digits) != ncol(df) | length(digits)==1) digits=rep(digits[1],ncol(df))
	ss=sapply(round(df),nchar)
	lar=df>1
	nn=sapply(names(df),nchar)
	
	# find longest string
	strw=sapply(1:ncol(df), function(x) max(nn, max(1,(ss[lar])+digits[x],na.rm=TRUE),na.rm=TRUE))  
	pr.df=data.frame(sapply(1:ncol(df), function(x) sprintf(paste("%",(strw[x]+buffer[x]),".",digits[x],"f",sep=""),df[,x])))   
	names(pr.df)=names(df)
	rownames(pr.df)=rownames(df)
	print(pr.df)
}
