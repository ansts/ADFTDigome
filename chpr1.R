#
# Read the densitometric data from microarrays. Needs a key.csv file 
# describing the files. bgcode is the code used for prescanned arrays (background)
#

chpr1=function(p, sh=FALSE, bgcode="ps", cntr=c("G","YPYDVPDYAG","KEVPALTAVETGAT"), chloc=NULL){
  require(parallel)
  require(future.apply)
  require(stringi)
  require(limma)
  require(reshape2)
  require(matrixStats)
  
  source("bgmy.R")
  
  lf=list.files(path = p)
  lf=lf[lf!="key.csv"]
  key=read.csv(paste(p,"key.csv", sep=""))
  N=nrow(key)
  
  # Get the files ----
  fls=lapply(lf,function(f){
    finp=paste(p,f,sep = "")
    fn=f
    fcon=file(finp)
    f2l=readLines(con=fcon, n=2)
    nlns=as.double(stri_extract(f2l[2], regex="[1-9]+"))
    fhead=readLines(con = fcon, n=nlns+2)
    x=read.delim(finp, skip=nlns+2, header=T, check.names = F, stringsAsFactors = FALSE)
    list(fn,x)
  })
  closeAllConnections()
  flns=lapply(fls,function(f){f[[1]]})
  fls=lapply(fls,function(f){f[[2]]})
  names(fls)=flns
  
  # Rearrange and summarize data ----
  bl=lapply(seq_along(key$file), function(i){
    f=fls[[key$file[i]]] 
    if (is.null(chloc)) {
      if (key$ch[i]=="G") 
          f=f[f$Block==key$block[i],c("Column","Row","ID","F532 Median","Flags")] 
      else if (key$ch[i]=="R") { 
          f=f[f$Block==key$block[i],c("Column","Row","ID","F635 Median","Flags")]
          colnames(f)=c("C", "R","ID","V","Fl")
      } else stop("Invalid detection channel identifier! Should be either G or R")
    }  else if (chloc=="G") 
          f=f[f$Block==key$block[i],c("Column","Row","ID","F532 Median","Flags")] 
      else 
          f=f[f$Block==key$block[i],c("Column","Row","ID","F635 Median","Flags")]
          colnames(f)=c("C", "R","ID","V","Fl")
     
    if (is.null(chloc)) chloc=key$ch[i]                                    
    return(list(key$pat[i],key$diag[i],chloc,f,key$file[i],key$PrePost[i])) #1-patient/2-diag/3-channel/4-[1-Col/2-Row/3-ID/4-F635 Med/5-Flags]/5-file/6 - prescan or assay
    })

  coor=lapply(bl, function(x){
    list(File=unlist(x[5]), Pat=x[1], Diagnosis=unlist(x[2]), Channel=unlist(x[3]), Row=x[4][[1]][,2], Column=x[4][[1]][,1])
    })                                                           #1-file/2-patient/3-diag/4-channel/5-row/6-column
  if (all(unlist(lapply(coor, function(x){unlist(coor[[1]][5:6])==unlist(x[5:6])})))) print("Ordered, homogeneous Sets") else return("Error: Nonhomogeneous Sets!")
  fl=lapply(bl,function(i){i[[4]][,5]})
  ID0=lapply(bl,function(i){i[[4]][,3]})
  #print("!")
  if( all(rowAlls(sapply(ID0, function(x) ID0[[1]]==x)))) ID0=ID0[[1]]

  pt=lapply(bl,function(i){as.character(i[[1]])})   #as.character(i[[3]])
  dg=lapply(bl,function(i){as.character(i[[2]])})
  ch=lapply(bl,function(i) i[[3]])
  PrePost=lapply(bl,function(i) i[[6]])
  cn=paste(PrePost,pt,dg, sep="_")
  
  Res0=as.data.frame(lapply(bl, function(x){x[[4]]$V}), col.names=cn)
  
  ij=aggregate(seq_along(fl), by=list(unlist(pt),unlist(dg)), function(x) as.numeric(x))
  ij=ij[order(ij$x[,1]),]
  R0=apply(ij,1,function(j){
    i=as.numeric(j[3:4])
    x=array(i, dimnames = list(unlist(PrePost)[i]))
    x=Res0[,x[names(x)=="post"]]-Res0[,x[names(x)=="pre"]]
    return(x)
  })
  cnm=apply(ij,1,function(j){
    paste(j[2],j[1], sep="_")
  })
  colnames(R0)=cnm
  # Reusing the name fls
  fls=unlist(lapply(coor,function(x){x$File}))
  flsR0=fls[ij$x[,1]]
 
  # R0=apply(R0,2,function(x){x=x-min(x)+1})          
  #rownames(R0)=ID
  #plan ("multisession", workers=16)
  B0=apply(R0, 2, function(clm){
      IM=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=clm)
      ID=data.frame(R=coor[[1]]$Row, C=coor[[1]]$Column, V=ID0, stringsAsFactors = FALSE)
      ID[,1]=as.double(ID[,1])
      ID[,2]=as.double(ID[,2])
      x=bgmy(IM,ID)
      x=x[order(x$R,x$C),]
      y=x[,1:2]
      x=x[,3]
      return(list(x,y))
  })                                                           
  closeAllConnections()
  
  B0xy=B0[[1]][[2]]
  B1=sapply(B0, function(l){l[[1]]})
  res=sapply(seq_along(R0[1,]),function(i){
    S=R0[,i]; B=B1[,i]
    Sr=S
    abZ=lm(Sr~B)
    j=cut(B,20,labels=F)
    bi=aggregate(B, by=list(j), "mean")
    z=abZ$residuals
    sdZ=sapply(sort(unique(j)), function(i) sd(z[j==i]))
    zz=lm(log10(sdZ)~log10(bi$x))
    a=zz$coefficients[2]
    b=zz$coefficients[1]
    sdf=(10^b)*(B^a)
    Sr=abZ$coefficients[2]*B+abZ$coefficients[1]+mean(abs(z))*z/sdf
    S=Sr
  res=backgroundCorrect.matrix(R0,B1, method = "normexp", normexp.method = "mle")
  
  return(list(coor,ID0,res,fl))
}

      plotarray <- function(X, Y, Z, nm) {
        cpl=colorRampPalette(c("#001030","#00507A","#007050","#9AAA00","#D0B000","#FF8000"))
        IM=data.frame(R=max(X)-X, C=Y, V=Z)
        y=acast(data=IM, C~R)
        par(mfrow=c(2,1))
        image(y, col=cpl(100), main=nm)
        par(mfrow=c(1,1))
        return(NULL)
      }