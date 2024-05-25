plotG3D=function(G, att="C", atfl=T, e_att=NULL, extatt=NULL, omit=NULL, 
                 coff=NULL, seed=777777, normaLap=T,smoothcol=T, interact=T){
  require(igraph)
  require(rgl)
  require(uwot)
  require(pbapply)
  
  cpl=colorRampPalette(c("#000000","#0050AA9F","#10AA109F","#FFFF009F","#FFA0009F","#B50000"), alpha=T)
  cpl4=c("#B5000090","#80808090","#10AA1090","#0020FF90")
  cm=components(G)
  G=induced.subgraph(G, V(G)[cm$membership==as.numeric(names(table(cm$membership)))[which.max(table(cm$membership))]])
  n=vcount(G)
  if (atfl) Att=vertex.attributes(G)[[att]] else Att=att[names(V(G))]
  if (normaLap) L=embed_laplacian_matrix(G,no=n-1, type = "I-DAD") else L=embed_laplacian_matrix(G,no=n-1)
  if( interact) {
    plot((n-50) : (n-1),L$D[(n-50) : (n-1)], xlab="# eigenvalue", ylab="The smallest non-zero eigenvalues")
    readline(" Cutoff OK?")
  }
  if (is.null(coff)) {
    dj=abs(diff(L$D))
    dj[(n-3):(n-2)]=0
    dj[1:(n/2)]=0
    j=(which.max(dj)+1):(n-1)
  } else j=(n-coff):(n-1)
  if (!is.null(omit[1]) & any(omit %in% j)) j=j[!(j %in% omit)] 
  print(j)
  Mx=L$X[,j]
  
  uMx=umap(Mx, n_neighbors = 50, min_dist = 0.005, spread=10,  n_epochs = 1500,
           verbose=F, n_components = 3, init = "normlaplacian",n_sgd_threads=1, seed=seed)
  
  #uMx=cmdscale(dist(Mx), k=7, add=T, list. = T, eig=T)
  # GOF=uMx$GOF
  # uMx=uMx$points

  if (smoothcol) vclr=cpl(vcount(G))[cut(Att, vcount(G), labels=F)] else vclr=Att
  names(vclr)=names(V(G))
  ecl=rep(1,ecount(G))
  if (!is.null(extatt)) eclr=extatt else {
      if (is.null(e_att)) {
        ee=ends(G, E(G))
        eclr=pbapply(ee,1,function(co) clrAvg(vclr[co[1]],vclr[co[2]]))
      } else {
        ecl=(edge_attr(G)[[e_att]]>4)*edge_attr(G)[[e_att]]
        eclr=ecl
        eclr[ecl==0|ecl>20]="#FFFFFF00"
        eclr[ecl>0&ecl<20]=cpl(20)[ecl]
      }
  }
    
    #eclr=cpl(ecount(G))[cut(edge_attr(G)[[e_att]], ecount(G), labels=F)]
 rglplot(G, layout=uMx, vertex.size=3, vertex.color=vclr, edge.color=eclr,
         edge.width=0.2*(ecl>0), vertex.label=NA, vertex.frame.color=NA)
 if (interact) readline("Adjust the window before hitting enter.")

  # pairs(Mx, col=vclr, pch=16, cex=0.5)
  return(list(zm=par3d()$zoom, usM=par3d()$userMatrix, wR=par3d()$windowRect))
}

clrAvg=function(co1,co2){
  require(broman)
  R1=hex2dec(substr(co1,2,3));G1=hex2dec(substr(co1,4,5));B1=hex2dec(substr(co1,6,7))
  R2=hex2dec(substr(co2,2,3));G2=hex2dec(substr(co2,4,5));B2=hex2dec(substr(co2,6,7))
  R3=dec2hex(floor((sqrt(R1)+sqrt(R2))/2)^2)
  G3=dec2hex(floor((sqrt(G1)+sqrt(G2))/2)^2)
  B3=dec2hex(floor((sqrt(B1)+sqrt(B2))/2)^2)
  if (nchar(R3)==1) R3=paste(0,R3,sep="", collapse="")
  if (nchar(G3)==1) G3=paste(0,G3,sep="", collapse="")
  if (nchar(B3)==1) B3=paste(0,B3,sep="", collapse="")
  return(paste("#",R3,G3,B3, collapse="", sep=""))
}