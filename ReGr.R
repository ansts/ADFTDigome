# Function building reactivity graph from a data matrix D . thr is the threshold
# of the criterion for establishing an edge

ReGr=function(D, thr=0.08){
  require(igraph)
  require(reshape2)
  require(parallel)
  require(pbapply)
  require(future.apply)
  require(stringdist)
  require(matrixStats)
  require(stringi)

  # Adjacency Matrix --------------------------------------------------------
  nc=14
  pep=rownames(D)
  Ds=D-min(D)+abs(min(D))/100                         
  ij=t(combn(length(pep),2))
  
  J=cut(seq_along(ij[,1]), nrow(ij) %/% 100, labels = F)
  plan("multisession", workers=nc)
  Fs=future_sapply(1:max(J), function(i){
    ii=which(J==i)   
    v=apply(ij[ii,],1,function(j){
      x=valD(Ds[j,])                                      
      return(x)
    })
    return(v)
  })
  closeAllConnections()
  Fs=unlist(Fs, recursive = F)
  # hist(Fs,breaks=1000)
  # lines(c(thr,thr),c(0,40), col=2)
  ij=rbind(ij,cbind(ij[,2],ij[,1]))
  x=data.frame(X=ij[,1], Y=ij[,2], V=c(Fs,Fs))
  ajmFs=acast(x,X~Y)
  rownames(ajmFs)=pep
  colnames(ajmFs)=pep

  #  Graph ------------------------------------------------------------------
  
  G=graph_from_adjacency_matrix(ajmFs*(ajmFs>thr), weighted = T, diag = F, mode = "undirected")
  G=simplify(G)
  return(G)
}

valF=function(m){
  require(reshape2)
  
  x=melt(m)
  2*var(aggregate(x$value, by=list(x$Var2), "mean")[,2])/mean(aggregate(x$value, by=list(x$Var2), "var")[,2])
  
}

vCCV=function(m){
  
  Cor=cor(t(m))[1,2]+1
  CV=sd(m)/mean(m)
  CVf=(1-CV/(0.3+CV))
  x=Cor+CVf    
  return(x)
  
}

valD=function(m){
  
  x=apply(m,1,function(l) quantile(l,c(0.25,0.5,0.75)))
  m=sapply(1:2, function(i) 2*(m[i,]-x[2,i])/(diff(x[c(1,3),i])))
  d=abs(rowDiffs(m))
  if (min(d)==0) d[d==min(d)]=min(d[d>0])/2
  x=quantile(d, c(0.2,0.8))
  d=d[d>x[1]&d<x[2]]
  res=(1/sum(d))
  return(res)
}