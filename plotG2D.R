plotG2D=function(G, att="K", e_att=NULL, omit=NULL, coff=NULL, seed=888888, normaLap=T, adj=F){
  require(igraph)
  require(rgl)
  require(uwot)
  
  cpl=colorRampPalette(c("#000000","#0050AA9F","#10AA109F","#FFFF009F","#FFA0009F","#B50000"), alpha=T)
  cm=components(G)
  G=induced.subgraph(G, V(G)[cm$membership==as.numeric(names(table(cm$membership)))[which.max(table(cm$membership))]])
  n=vcount(G)
  if (normaLap) L=embed_laplacian_matrix(G,no=n-1, type = "I-DAD") else L=embed_laplacian_matrix(G,no=n-1)
  if(adj) {
    plot(L$D[(n-50) : (n-1)])
    readline(" Cutoff OK?")
  }
  if (is.null(coff)) {
      dj=abs(diff(L$D))
      dj[(n-3):(n-2)]=0
      dj[1:(n/2)]=0
      j=(which.max(dj)+1):(n-1)
  } else j=(n-coff):(n-1)
  if (!is.null(omit[1]) & any(omit %in% j)) j=j[!(j %in% omit)] 
  Mx=L$X[,j]
  uMx=umap(Mx, n_neighbors = 50, min_dist = 0.005, spread=10,  n_epochs = 1500, 
           verbose=F, n_components = 2, init = "normlaplacian",n_sgd_threads=1, seed=seed)
  vclr=cpl(vcount(G))[rank(vertex.attributes(G)[[att]])]
  names(vclr)=names(V(G))
  if (is.null(e_att)) {
    ee=ends(G, E(G))
    eclr=pbapply(ee,1,function(co) clrAvg(vclr[co[1]],vclr[co[2]]))
  } else eclr=cpl(ecount(G))[cut(edge_attr(G)[[e_att]], ecount(G), labels=F)]
  plot(G, layout=uMx, vertex.size=3, vertex.color=vclr, edge.color=eclr, edge.width=0.15, vertex.label=NA)
  return(j)
}