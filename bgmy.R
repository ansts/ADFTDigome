bgmy=function(IM,ID){
  require(reshape2)
  require(e1071)
  
  cpl=colorRampPalette(c("#001030","#0050AA","#10AA10","#FFFF00","#FFA000","#B50000"))
  img=acast(IM, names(IM)[1]~names(IM)[2], value.var = names(IM)[3])
  
  nr=nrow(img)
  nc=ncol(img)
  id=unique(ID)
  
  # Find targets that are represented by more than two spots,
  # create all couples of them with surrogate IDs so that strictly couples of 
  # spots can be used further but at the same time all available couples are used
  
  tID=table(ID[,3])
  nmrs=names(tID[tID>2])
  if (length(nmrs)>0){
    nm=length(nmrs)
    nmr=lapply(1:nm,function(i){cbind(ID[ID[,3]==nmrs[i],],IM[ID[,3]==nmrs[i],3])})
    nmr=lapply(nmr,function(i){colnames(i)=c("R","C","ID","V"); i})
    nmri=unlist(lapply(nmr,function(n){which(ID[,3] %in% n[,3])}))
    ID=ID[-nmri,]
    IM=IM[-nmri,]

    nmr=lapply(nmr,function(n){n[combn(nrow(n),2),]})
    nmr=lapply(nmr,function(n){n[,3]=paste(n[,3],"_",rep(1:(nrow(n) %/% 2),each=2),sep = ""); n})
    for (i in 1:nm) {
      x=nmr[[i]][,1:3]
      colnames(x)=colnames(ID)
      ID=rbind(ID,x)
    }
    for (i in 1:nm) {
      x=nmr[[i]][,c(1,2,4)]
      colnames(x)=colnames(IM)    
      IM=rbind(IM,x)
    }
  }
  B=IM

  # The algorithm
  
  ds=lapply(seq_along(ID[,1]), function(rx){
    pn=ID[rx,3]
    rc=as.double(ID[rx,1:2])
    y=B[ID[,3]==pn,]
    p1=y[y[,1]==rc[1] & y[,2]==rc[2],]
    p2=y[y[,1]!=rc[1] | y[,2]!=rc[2],]
    dx=p1[,3]-p2[,3]
    ff=c(which(ID[,1]==p1[,1] & ID[,2]==p1[,2] & ID[,3]==pn),which(ID[,1]==p2[,1] & ID[,2]==p2[,2] & ID[,3]==pn)) # indices of the couple
    #print(ff)
    c(rc,dx,ff)
  })                                       #  differences
  ds=t(as.data.frame(ds))
  colnames(ds)=c("R","C","V","P1","P2")
  rownames(ds)=ID[,3]
  fl1=TRUE
  co=1
  bas=list()
  bas[1]=1
  
  repeat {
    if (fl1==TRUE) B[,3]=ds[,3]-mean(ds[,3])
    fl1=FALSE
 
    Av=t(sapply(seq_along(IM[,1]),function(x){
        r=B[x,1]
        c=B[x,2]
        rmax=max(B[,1])
        cmax=max(B[,2])
        
        tr=(r-1):(r+1)
        tc=(c-1):(c+1)
        if (r==1) tr=1:3
        if (c==1) tc=1:3
        if (r==rmax) tr=(rmax-3):rmax
        if (c==cmax) cr=(cmax-3):cmax
      
        dsi=B[B[,1] %in% tr & B[,2] %in% tc,]
        if (nrow(dsi)>9){
          dsi=t(apply(unique(dsi[,1:2]),1, function(x){y=mean(dsi[dsi[,1]==x[1] & dsi[,2]==x[2],3]); c(x,y)}))
          colnames(dsi)==c("R","C","V")
        }
        
        dsi=c(r,c,mean(dsi[,3]))   # mean difference in the  vicinity
        return(dsi)
      }))
    B[,3]=Av[,3]
    #lmd=svm(B[,1:2], B[,3], cost = 500, gamma=4)   
    #B[,3]=predict(lmd, B[,1:2])
   
    bs=list()
    B0=B
    
    for (i in seq_along(B[,1])){
      jr=ds[ds[i,5],1]
      jc=ds[ds[i,5],2]
      ix=B[,1]==jr & B[,2]==jc & ID[,3]==ID[i,3]
      bs[i]=B[ix,3]-B[i,3]+ds[i,3]                        # 
      B0[ix,3]=(5*B[ix,3]+B[i,3]-ds[i,3])/6             # Replace the partner value with ..
    }
    B=B0
    co=co+1
    bas[co]=mean(unlist(bs[bs>0]))
    crit=abs(bas[[co]]-bas[[co-1]])/bas[[co-1]]
    print(c(bas,crit))
    if (crit<0.0000001) break

  }
  lmd=svm(B[,1:2], B[,3], cost = 1000, gamma=0.4)   
  #lmd0=svm(IM[,1:2], IM[,3], cost = 1000, gamma=0.4)
  Bsm=predict(lmd, B[,1:2])
  Bsm=aggregate(Bsm, by=list(Row=B$R, Column=B$C), FUN=mean)
  colnames(Bsm)=colnames(IM)
  Bsm[,3]=Bsm[,3]-min(Bsm[,3])
  #Bsm=aggregate(B[,3], by=list(Row=B$R, Column=B$C), FUN=mean)
  return(Bsm)
}