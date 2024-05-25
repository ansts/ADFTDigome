adjL2G=function(AL){
require(igraph)
require(reshape2)

L0=lapply(AL, as.list)
L=L0[lengths(L0)>0]
if (length(L)<4) {
  print("Less than 4 verteces")
  return(NULL)
}
EL=melt(L)
EL=EL[,c(3,1,2)]
EL=t(apply(EL,1,function(r){
  r1=r[1:2]
  r1=r1[order(r1)]
  return(r1)
}))

EL=unique(EL)
rownames(EL)=seq_along(EL[,1])
colnames(EL)=c("From", "To")
EL=as.data.frame(EL,stringsAsFactors=F)

G=graph_from_data_frame(EL, directed=F)
G=simplify(G)
#plot(G, vertex.size=1, main=paste("Graph at Cut Off",coff, collapse = " "))
return(G)
}