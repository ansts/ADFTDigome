# Load packages ------
require(limma)
require(dunn.test)
require(pbapply)
require(prallel)
require(stringdist)
require(reshape2)
require(matrixStats)
require(igraph)
require(chisq.posthoc.test)
require(gplots)
require(vioplot)
require(future.apply)
require(aricode)
require(stringi)
require(factoextra)
require(FactoMineR)
require(Rfast)
require(uwot)
require(corrplot)
require(rBLAST)
require(patchwork)
require(rgl)
require(pandoc)
require(htmltools)
require(htmlwidgets)
require(ppcor)
require(biomartr)
require(easyalluvial)

fxs=paste(getwd(),c("chpr1.R", "bgmy.R","mkZpepMx.R", "pepnorm.R", "ReGr.R","plotG2D.R",
                    "plotG3D.R","rfe.R","clucri.R"),
          sep="/")
sapply(fxs,source)

# Some databases ------

load("IgJT")     # human Ig J regions from NCBI 
Ig7=t(qgrams(IgJtrim, q=7))
Ig7F=Ig7[,1]
names(Ig7F)=rownames(Ig7)
Ig7=rownames(Ig7)
idnnL=pbsapply(pp, function(p){
  Ig7F[stringdist(p,Ig7, method = "lcs")<3] 
})
idnn=lengths(idnnL)
idnnt=table(idnn)

EBVepi=read.csv("EBVautoImmEp..csv")
EBVepi=unlist(stri_extract_all(EBVepi[,1], regex="^\\w+"))
EBVepi7=qgrams(EBVepi,q=7)[1,]

load("matpep")  # 7-mer Ph.D.-7 sequences obtained after single amplification without selection 
load("ODxs790_30.out")  # healthy donor public IgM mimotope library
orilab=unlist(lapply(ODxs790_30.out, function(Tb){
  Tb$Alignment
}))


# Analysis of the IgG and IgM Binding ------------------------------------------

cpl=colorRampPalette(c("#000000","#0050AA9F","#10AA109F","#FFFF009F","#FFA0009F","#B50000"), alpha=T)


## Collect the data from the .gpr files in a folder /gpr 
#  containing also a key.csv file describing the slides 

pth="bigarraydata//"
dgn=rep(rep(c("A","F","ND","C"),5),2)
pat=rep(c("1.6","1.6","1.6","1.3","2.7","2.7","2.7","2.6","3.8","3.8","3.8","3.4","4.8","4.1","4.1","4.2","5.2","5.2","5.2","5.1"),2)
lf=list.files(path = pth)
PrePost=c(rep("post",20),rep("pre",20))
fnms=lf
key=data.frame(file=fnms,pat=pat,block=1,diag=dgn,PrePost=PrePost,ch="both")
write.csv(key, file=paste(pth,"key.csv", sep=""))

blgr=read.csv(file="blGr1.csv")
blgr=blgr[,-1]
colnames(blgr)=c("C","A","F","ND")

# chpr1 extracts the data from the gpr files and organizes it subtracting the
# background. A list of two matrices is created - the first with the IgM data,
# the second - with the IgG data 
nms=c("IgM","IgG")

DnGM=lapply(nms, function(ch){
  if (ch=="IgG") ch="R"
  if (ch=="IgM") ch="G"
  alldata=chpr1(pth, sh=T, chloc=ch)
  Coor=alldata[[1]] # file/patients/diagnosis/channel/row/column
  iD=alldata[[2]]   # ID of spots
  Res=alldata[[3]]  # Results - locally normalized data (the background is reconstituted using the duplicate diff. and subtracted using normexp)
  rownames(Res)=iD  # 
  Fl=alldata[[4]]   # Flags data
  
  FLdf=as.data.frame(Fl[1:ncol(Res)])
  FL=rowSums(FLdf)==0
  
  # WD - data filtered for flags from the densitormetry then aggregated to average
  # the duplicate peptide spots
  WD=Res[FL,]
  WDa=aggregate.data.frame(WD, by=list(rownames(WD)), FUN=mean)
  WDam=data.frame(log10(WDa[,2:ncol(WDa)]), stringsAsFactors = FALSE, row.names = WDa[,1])
  
  #### normalize for amin oacid composition dependent binding 
  #(e.g. non-specific stickiness of charged amino acids)
  
  D=pepnorm(WDam)
  
  #### normalize between arrays within topology/isotype groups 
  Dn=normalizeCyclicLoess(D, method="affy", iterations=3) 
  return(Dn)
})  

names(DnGM)=nms
ppM=rownames(DnGM$IgM)
lM=length(ppM)
ppG=rownames(DnGM$IgG)
lG=length(ppG)
setdiff(ppM,ppG)
setdiff(ppG,ppM)
pp=ppM; rm(ppG,ppM)
diapat=t(sapply(colnames(DnGM$IgM), function(y) unlist(strsplit(y, split="_|\\."))))
diapat=data.frame(dia=diapat[,1], pat1=as.numeric(diapat[,2]), pat2=as.numeric(diapat[,3]))
diapat=cbind(1:20,diapat)
dia=diapat$dia
jup=lapply(seq_along(diapat[,1]),function(i){
  l=diapat[i,]
  j=which(diapat[,"dia"]==l[[2]]&(diapat[,"pat2"]==l[[3]]|diapat[,"pat2"]==l[[4]]))
  return(j)
}) 
names(jup)=rownames(diapat)

# Blood group variables ------

zblgr=t(apply(diapat,1,function(l) {
  x=paste(blgr[as.numeric(unlist(l[3:4])),l[[2]]], sep = "", collapse = "")
  x=strsplit(x, split="")[[1]]
  x=x[x %in% c("A","B")]
  t0=rep(0,2)
  names(t0)=c("A","B")
  if (length(x)>0) {
    t1=table(x)
    t0[names(t1)]=t1
  }
  return(t0)
}))
diapat=cbind(diapat,zblgr)
diall=as.numeric(as.factor(diapat$dia))

blAB=lapply(DnGM, function(M){   # correlation of each peptide reactivity with A/B expression
  t(apply(M,1,function(l) {
    apply(diapat[,5:6], 2, function(clm){
      lmx=summary(lm(l~clm))
      return(lmx$coefficients[2,1])
    })
  }))
})

blABp=lapply(DnGM, function(M){
  t(apply(M,1,function(l) {
    p=apply(diapat[,5:6], 2, function(clm){
      lmx=summary(lm(l~clm))
      return(lmx$coefficients[2,4])
    })
    p.adjust(p)
  }))
})

## Fig  -------

rM=quantile(DnGM$IgM,c(0.003,0.997))
rG=quantile(DnGM$IgG,c(0.003,0.997))

pdf(file="rawdata.pdf", width=10,height=10)
  pairs(DnGM$IgM, cex=0.05,col=rgb(0,0,0,0.2), main="IgM", xlim=rM, ylim=rM)
  pairs(DnGM$IgG, cex=0.05,col=rgb(0,0,0,0.2), main="IgG", xlim=rG, ylim=rG)
dev.off()

# Correlations between diagnosis means -----

diaMeans=lapply(nms, function(iso){
  mns=rowmeans(DnGM[[iso]])
  X=lm(DnGM[[iso]]~rowMeans(DnGM[[iso]]))
  t(apply(X$residuals,1,function(l){
    x=aggregate(l, by=list(diapat$dia), "mean")
    n=x$Group.1; x=x$x; names(x)=n
    return(x)
  }))
})

# diaMeans=lapply(nms, function(iso){
#   t(apply(DnGM[[iso]],1,function(l){
#     x=aggregate(l, by=list(diapat$dia), "mean")
#     n=x$Group.1; x=x$x; names(x)=n
#     return(x)
#   }))
# })

names(diaMeans)=nms

cpl2=colorRampPalette(c("#000000","#0050AA9F","#FFFFFF","#FFA0009F","#B50000"), alpha=T)
diaMnside=cbind(diaMeans[[1]][pp,], diaMeans[[2]][pp,],idnn[pp], blAB$IgM[pp,1],blAB$IgM[pp,2],blAB$IgG[pp,1],blAB$IgG[pp,2])
colnames(diaMnside)=c(paste(colnames(diaMnside)[1:8],c(rep("IgM",4),rep("IgG",4)), sep=":"),"Id NN", 
                      "BlGr A/M", "BlGr B/M","BlGr A/G", "BlGr B/G")
corrplot(pcor(diaMnside)$estimate, method="color", order="hclust",hclust.method = "ward.D2", col=cpl2(15), diag=F)
corrplot(pcor(diaMnside)$estimate, method="color", order="hclust",hclust.method = "ward.D2", col=cpl2(15), diag=F)

##### Fig -----
X=cor(diaMnside)
# X=X$estimate
# colnames(X)=colnames(diaMnside)
# rownames(X)=colnames(diaMnside)
corrplot(X, method="color", order="hclust",hclust.method = "ward.D2", col=cpl2(15), diag=F)
heatmap.2(X, hclustfun = function(d) hclust(d, method="ward.D2"), col=cpl2(15), trace="none")

# PCA -------------------------------------------------------------------

pcaDn=lapply(DnGM, function(M){
  PCA(M,ncp=20)
})
mnDn=lapply(DnGM, rowMeans)
apply(pcaDn$IgM$var$cor, 2, function(cl) {
  x=aggregate(cl, by=list(diapat$dia), "mean")
  nx=x$Group.1
  x=x$x;names(x)=nx
  return(x)
})
cor(mnDn$IgM,pcaDn$IgM$ind$coord)
cor(mnDn$IgG,pcaDn$IgG$ind$coord)


cpl1=c("red","black","green","blue")


im=c(1,2,4,8); iim=combn(im,2)[,c(1,4,5,6)]
ig=c(1,2,3,5,10); iig=combn(ig,2)[,c(1,5,6,7)]
ii=list(iim,iig)

dimmaxM=rowMaxs(diaMeans$IgM)
dimmaxG=rowMaxs(diaMeans$IgG)
dimmax=list(dimmaxM, dimmaxG)

fvizplots=lapply(1:2, function(i){
  apply(ii[[i]],2,function(ax) {
    fviz_pca_var(pcaDn[[i]], 
                 axes=ax, palette=cpl1, col.var = as.factor(diapat$dia), label="none")})
})
(fvizplots[[1]][[1]]+fvizplots[[1]][[2]])/(fvizplots[[1]][[3]]+fvizplots[[1]][[4]])
(fvizplots[[2]][[1]]+fvizplots[[2]][[2]])/(fvizplots[[2]][[3]]+fvizplots[[2]][[4]])

dimmax=list(dimmaxM, dimmaxG)

fvizplotsx=lapply(1:2, function(i){
  apply(ii[[i]],2,function(ax) {
    fviz_pca_ind(pcaDn[[i]], 
                 axes=ax, palette=cpl1, 
                 col.ind = as.factor(c("A","C","F","ND")[dimmax[[i]]]), 
                 label="none")})
})

(fvizplots[[1]][[1]]+fvizplotsx[[1]][[2]])/(fvizplotsx[[1]][[3]]+fvizplotsx[[1]][[4]])
(fvizplots[[2]][[1]]+fvizplotsx[[2]][[2]])/(fvizplotsx[[2]][[3]]+fvizplotsx[[2]][[4]])

corpcaM=cor(diaMeans$IgM,pcaDn$IgM$ind$coord[rownames(diaMeans$IgM),])
heatmap.2(corpcaM[,c(2:8)], hclustfun = function(d) hclust(d,method = "complete"), col=cpl(100), trace="none")
corpcaG=cor(diaMeans$IgG,pcaDn$IgG$ind$coord[rownames(diaMeans$IgG),])
heatmap.2(corpcaG[,c(2:8)], hclustfun = function(d) hclust(d,method = "complete"), col=cpl(100), trace="none")



##### Fig  ------------
 
(fvizplots[[1]][[1]]+fvizplotsx[[1]][[2]])/(fvizplots[[2]][[1]]+fvizplotsx[[2]][[2]])

M=data.frame(IgM=(DnGM$IgM-min(DnGM$IgM)/diff(range(DnGM$IgM))), 
             IgG=(DnGM$IgG-min(DnGM$IgG)/diff(range(DnGM$IgG))))



####


pcaDnGM=PCA(M,ncp=20)

fviz_screeplot(pcaDnGM, ncp=20)

jj=list(c(1,2),c(2,3),c(3,4),c(4,5))
jointpca=lapply(1:4,function(j){
    fviz_pca_var(pcaDnGM, axes=jj[[j]], repel=T, 
                 col.var=paste(rep(diapat$dia,2),rep(c("M","G"), each=20),sep="|"), 
                 col.label=paste(rep(diapat$dia,2),rep(c("M","G"), each=20),sep="|"),
                 labelsize=0.4)
})

(jointpca[[1]]+jointpca[[2]])/(jointpca[[3]]+jointpca[[4]])

corjoint=cor(diaMnside,pcaDnGM$ind$coord[rownames(diaMnside),])
heatmap.2(corjoint[,-1], hclustfun = function(d) hclust(d,method = "ward.D2"), col=cpl(100), trace=NULL)

# ReaGraph ---------------------------------------------------------------------

mn=sapply(DnGM,rowMeans)
std=sapply(DnGM, rowSds)
pWRST=lapply(DnGM,function(iso){
  t(apply(iso,1,function(ro){
    z=aggregate(rank(ro), by=list(diapat$dia), "sum")
    x=z$x; names(x)=z$Group.1; z=x
    dun=dunn.test(ro, diapat$dia, method="none")
    Z=dun$Z; names(Z)=dun$comparisons
    P=dun$P.adjusted; names(P)=paste(names(Z),"_padj",sep="_")
    return(c(z,Z,P))
  }))
})
X=pWRST$IgM[,11:16]
X=apply(X,2,p.adjust, method="BH")
pWRST$IgM[,11:16]=X
X=pWRST$IgG[,11:16]
X=apply(X,2,p.adjust, method="BH")
pWRST$IgG[,11:16]=X

pMndia=lapply(DnGM,function(iso){
  t(apply(scale(iso),1,function(ro){
    z=aggregate(ro, by=list(diapat$dia), "mean")
    x=z$x
    names(x)=z$Group.1
    return(x)
  }))
})

reabydia=lapply(nms, function(iso){
  t(apply(t(scale(t(DnGM[[iso]]))),1,function(l){
    x=aggregate(l, by=list(diapat$dia), "mean")
    n=x$Group.1; x=x$x; names(x)=n
    return(x)
  }))
})
names(reabydia)=nms

# Scrambled matrices for bootstrapping different statistics

set.seed(13071962)

DnGMBS=lapply(nms, function(iso){
  print(iso)
  pblapply(1:100, function(i){
    array(sample(DnGM[[iso]]), dim=dim(DnGM$IgM), dimnames = dimnames(DnGM[[iso]]))
  })
})
names(DnGMBS)=nms

reabydiaBS=lapply(nms, function(iso){
  pblapply(1:100, function(z){
    t(apply(t(scale(t(DnGMBS[[iso]][[z]]))),1,function(l){
      x=aggregate(l, by=list(diapat$dia), "mean")
      n=x$Group.1; x=x$x; names(x)=n
      return(x)
    }))
  })
})
names(reabydiaBS)=nms

### Graphs by Isotype ----------------------

# the actual graphs are in the list Gi
Gi=lapply(nms, function(iso){
  g=ReGr(DnGM[[iso]])
  g=set_vertex_attr(g, name="Mean", value=mn[,iso])
  g=set_vertex_attr(g, name="SD", value=std[,iso])
  g=set_graph_attr(g,name="name", value=iso)
  
  xm=pWRST[[iso]][,1:10]
  for (i in 1:10) g=set_vertex_attr(g, name=colnames(xm)[i], value=xm[,i])
  return(g)
})
names(Gi)=sapply(Gi, function(g) graph_attr(g)$name)

Gi[[1]]=set_vertex_attr(Gi[[1]], name="blgrA", value=blAB[[1]][,1])
Gi[[1]]=set_vertex_attr(Gi[[1]], name="blgrB", value=blAB[[1]][,2])
Gi[[2]]=set_vertex_attr(Gi[[2]], name="blgrA", value=blAB[[2]][,1])
Gi[[2]]=set_vertex_attr(Gi[[2]], name="blgrB", value=blAB[[2]][,2])
Gi[[1]]=set_vertex_attr(Gi[[1]], name="blgrAsig", value=(blABp[[1]][,1]<0.05)*sign(blAB[[1]][,1]))
Gi[[1]]=set_vertex_attr(Gi[[1]], name="blgrBsig", value=(blABp[[1]][,2]<0.05)*sign(blAB[[1]][,2]))
Gi[[2]]=set_vertex_attr(Gi[[2]], name="blgrAsig", value=(blABp[[2]][,1]<0.05)*sign(blAB[[2]][,1]))
Gi[[2]]=set_vertex_attr(Gi[[2]], name="blgrBsig", value=(blABp[[2]][,2]<0.05)*sign(blAB[[2]][,2]))
Gi$IgM=set_vertex_attr(Gi$IgM, name="idnn", value=log10(idnn[names(V(Gi$IgM))]+1))
Gi$IgG=set_vertex_attr(Gi$IgG, name="idnn", value=log10(idnn[names(V(Gi$IgG))]+1))

## Exchange labels -----------------------

Gi$IgM=set_vertex_attr(Gi$IgM, name="A_g", value=vertex_attr(Gi$IgG)$A)
Gi$IgM=set_vertex_attr(Gi$IgM, name="F_g", value=vertex_attr(Gi$IgG)$F)
Gi$IgM=set_vertex_attr(Gi$IgM, name="ND_g", value=vertex_attr(Gi$IgG)$ND)
Gi$IgM=set_vertex_attr(Gi$IgM, name="C_g", value=vertex_attr(Gi$IgG)$C)
Gi$IgG=set_vertex_attr(Gi$IgG, name="A_m", value=vertex_attr(Gi$IgM)$A)
Gi$IgG=set_vertex_attr(Gi$IgG, name="F_m", value=vertex_attr(Gi$IgM)$F)
Gi$IgG=set_vertex_attr(Gi$IgG, name="ND_m", value=vertex_attr(Gi$IgM)$ND)
Gi$IgG=set_vertex_attr(Gi$IgG, name="C_m", value=vertex_attr(Gi$IgM)$C)

plot(degree(Gi$IgM),degree(Gi$IgG), pch=16, col=rgb(0,0,0,0.5))

sapply(Gi, components)

maxclq=pblapply(Gi, function(g){
  max_cliques(g)
})

#### GM crossreactivities --------

GMCR=sapply(pp, function(p){
  lm=DnGM$IgM[p,]
  lg=DnGM$IgG[p,]
  valD(rbind(lm,lg))
})

Gi$IgM=set_vertex_attr(Gi$IgM, name="GMCR", value=log10(GMCR[names(V(Gi$IgM))]))
Gi$IgG=set_vertex_attr(Gi$IgG, name="GMCR", value=log10(GMCR[names(V(Gi$IgG))]))

assortativity(Gi$IgM,vertex_attr(Gi$IgM)$GMCR)
assortativity(Gi$IgG,vertex_attr(Gi$IgG)$GMCR)

save(Gi, file="Gi")

## Singletons  -----

sngltn=lapply(nms, function(iso) {
  x=which(components(Gi[[iso]])$membership!=1)
  jn=names(x)
  m=as.data.frame(vertex_attr(Gi[[iso]])[4:7])[names(V(Gi[[iso]])) %in% jn,]
  rownames(m)=jn
  heatmap.2(as.matrix(m), dendrogram = "both", col=cpl(100), 
            hclustfun = hclust, method="ward.D2", trace="none")
  vioplot(m, names=colnames(m), h=4)
  return(m)
})

## Bootstrap graph ------

GiBS=lapply(nms, function(iso){
  pblapply(DnGMBS[[iso]], function(M){
    g=ReGr(M)
    pWRST=t(apply(M,1,function(ro){
      z=aggregate(rank(ro), by=list(diapat$dia), "sum")
      x=z$x; names(x)=z$Group.1; z=x
      dun=dunn.test(ro, diapat$dia, method="none")
      Z=dun$Z; names(Z)=dun$comparisons
      P=dun$P.adjusted; names(P)=paste(names(Z),"_padj",sep="_")
      return(c(z,Z,P))
    }))
    xm=pWRST[,1:10]
    for (i in 1:10) g=set_vertex_attr(g, name=colnames(xm)[i], value=xm[,i])
    if (iso=="IgM") { 
      g=set_vertex_attr(g, name="blgrA", value=blAB[[1]][,1])
      g=set_vertex_attr(g, name="blgrB", value=blAB[[1]][,2])
    } else {
      g=set_vertex_attr(g, name="blgrA", value=blAB[[2]][,1])
      g=set_vertex_attr(g, name="blgrB", value=blAB[[2]][,2])
    }
    return(g)
  })
})
names(GiBS)=nms
GiBS=lapply(nms, function(iso){
  isalt=nms[nms!=iso]
  lapply(1:100, function(i){
    g0=GiBS[[iso]][[i]];g1=GiBS[[isalt]][[i]]
    g0=set_vertex_attr(g0, name="A_a", value=vertex_attr(g1)$A)
    g0=set_vertex_attr(g0, name="F_a", value=vertex_attr(g1)$F)
    g0=set_vertex_attr(g0, name="ND_a", value=vertex_attr(g1)$ND)
    g0=set_vertex_attr(g0, name="C_a", value=vertex_attr(g1)$C)
    return(g0)
  })
})
names(GiBS)=nms
save(DnGMBS,GiBS, file="GiBSdata")

#### Components  -----

grcomp=sapply(Gi, function(g) max(components(g)$csize))

compBS=lapply(nms, function(iso){
  sapply(GiBS[[iso]], function(g) table(components(g)$csize))
})

compBSlcm=sapply(compBS, function(L){
  sapply(L, function(l) max(as.numeric(names(l))))
})
colnames(compBSlcm)=nms

zcomp=sapply(1:2, function(i){
  X=compBSlcm[,i]
  (max(components(Gi[[i]])$csize)-mean(X))/sd(X)
})

###### Density  ====

bsdns=sapply(1:2, function(i){
  sapply(GiBS[[i]], graph.density)
})

zdns=sapply(1:2, function(i){
  (graph.density(Gi[[i]])-mean(bsdns[,i]))/sd(bsdns[,i])
})

## Assortativities ----------------------------

assrt=lapply(c("IgM","IgG"), function(iso){
  attt=vertex_attr_names(Gi[[iso]])[c(4:7,14:15,18:23)]
  x=sapply(attt, function(x) {
    assortativity(Gi[[iso]], vertex_attr(Gi[[iso]])[[x]])
  })
  nch=sub("_\\w+", "_a", names(x))
  names(x)=nch
  return(x)
})
names(assrt)=c("Assortativity IgM","Assortativity IgG")

assGiBS=lapply(nms, function(iso){
  gbs=GiBS[[iso]]
    t(pbsapply(gbs, function(g){
      g=induced.subgraph(g, V(g)[names(V(Gi[[iso]]))])
      attt=vertex_attr_names(Gi[[iso]])[c(4:7,14:15,18:22,25)]
      x=sapply(attt, function(x) {
        assortativity(g, vertex_attr(Gi[[iso]])[[x]])
      })
      nch=sub("_\\w+", "_a", names(x))
      names(x)=nch
      return(x)
    }))
  })

names(assGiBS)=c("IgM","IgG")
par(mfrow=c(1,2))
for (i in 1:2){
  boxplot(assGiBS[[i]][,c(1:7,10,8,9,11,12)], ylim=c(-1,1), xlab="", ylab="Assortativity", las=2)
  par(new=T)
  plot(assrt[[i]][c(1:7,10,8,9,11,12)], xlim=c(0.5,12.5), ylim=c(-1,1), cex=1.5, pch=16, col=2, xaxt="n", yaxt="n",xlab="", ylab="", main=nms[i])
}


## Cliques ----

maxclqBS=lapply(GiBS, function(gi){
  pbsapply(gi, function(g){
    table(lengths(max_cliques(g)))
  })
})

# Idiotypy --------------------

idnnPHDLBS=lapply(1:100, function(i){
  print(i)
  pbsapply(sample(matpep,722), function(p){
    Ig7F[stringdist(p,Ig7, method = "lcs")<3] 
  })
})  
idnnPHDBSt=lapply(idnnPHDLBS, function(L) table(lengths(L)))
save(idnnPHDLBS, file="idnnPHDLBS")

idnnPHDBStn=names(table(unlist(sapply(idnnPHDBSt, function(x) as.numeric(names(x))))))

idnnIgMLBS=pblapply(1:100, function(i){
  print(i)
  sapply(sample(orilab,722), function(p){
    Ig7F[stringdist(p,Ig7, method = "lcs")<3] 
  })
}) 
idnnIgMBSt=lapply(idnnIgMLBS, function(L) table(lengths(L)))
save(idnnIgMLBS, file="idnnIgMLBS")

idnnIgMBStn=names(table(unlist(sapply(idnnIgMBSt, function(x) as.numeric(names(x))))))

grtt=rep(0,366)
names(grtt)=0:365
for (i in 1:100){
  grtt[names(idnnIgMBSt[[i]])]=grtt[names(idnnIgMBSt[[i]])]+idnnIgMBSt[[i]]
  grtt[names(idnnPHDBSt[[i]])]=grtt[names(idnnPHDBSt[[i]])]+idnnPHDBSt[[i]]
} 
grtt[names(idnnt)]=grtt[names(idnnt)]+idnnt
grttSum(grtt)
jx=c(0,1,2,3,4,6,9,14,31,365)
idnnIgMBStAgg=t(sapply(idnnIgMBSt,grttSum))
idnnPHDBStAgg=t(sapply(idnnPHDBSt,grttSum))
TTT=lapply(as.character(jx), function(n) cbind(idnnPHDBStAgg[,n],idnnIgMBStAgg[,n]))
x=TTT[[1]]
for (i in 2:10){
  x=cbind(x,TTT[[i]])
}
TTT=x
colnames(TTT)=paste(c("PHD", "IgM"),rep(jx,each=2), sep = "_")

boxplot(TTT, ylim=c(0,250), las=2, col=rep(c("#F0F0F0","#A0A0A0"),10), 
        main="Distribution of number of sequences with idiotope NN", 
        ylab="Number of sequences")
title(xlab = "Number of NN - upper limit of the interval", line=4.5)
par(mar=c(8, 4.1,4.1, 2.1))
points((1:10)*2-0.5,grttSum(idnnt), ylim=c(0,250), pch=16, cex=2, col=2)
legend("topright", 
       legend=c("Ph.D.-7 library", "IgM mimotopes", "Selected mimotopes"), 
       fill=c("#F0F0F0","#A0A0A0", "red"), bty = "n")


plotG3D(Gi$IgM, att="idnn", coff=4)
plotG3D(Gi$IgG, att="idnn", coff=4)

idassBS=sapply(nms, function(iso){
  pbsapply(1:1000, function(i){
    idx=sample(vertex_attr(Gi[[iso]])$idnn)
    assortativity(Gi[[iso]], idx)
  })
})

asidG=assortativity(Gi$IgG, idnn[names(V(Gi$IgG))])
asidM=assortativity(Gi$IgM, idnn[names(V(Gi$IgM))])

ecdfM=ecdf(idassBS[,1])
ecdfG=ecdf(idassBS[,2])
asidqM=ecdfM(asidM)
asidqG=ecdfG(asidG)

boxplot(idassBS, notch=T, outline=F, ylim=c(-0.04,0.04), ylab="Id assortativity")
points(1:2,c(asidM, asidG), col=2, pch=16, cex=2)
text(c(1.175,2.175), c(0.035,0.02), labels=c(asidqM,asidqG))

attidcor=t(sapply(nms, function(iso){
  Mx=cbind(as.data.frame(reabydia[[iso]][names(V(Gi[[iso]])),]),as.data.frame(vertex_attr(Gi[[iso]])[c(4:7,14:15,22)]))
  cor(vertex_attr(Gi[[iso]])$idnn,Mx)
}))
colnames(attidcor)=c(paste("mn",colnames(reabydia[[iso]]),sep="_"), vertex_attr_names(Gi$IgM)[c(4:7,14:15,22)])

xi=as.numeric(names(table(idnn)))
xi=xi[-length(xi)]
WRSTidnn=lapply(nms, function(iso){
  plot(NULL,xlim=c(1,100), ylim=c(-7,7), ylab="z")
  pbsapply(1:4, function(j){
    z=sapply(xi, function(i){
      fc=cut(idnn, c(min(idnn)-1,i,max(idnn)+1),labels=F)
      Y1=reabydia[[iso]][fc==2,j]
      Y2=reabydia[[iso]][fc==1,j]
      res=wilcox.test(x=Y1,y=Y2)
      R=res$statistic
      n1=length(Y1)
      n2=length(Y2)
      mu=(n1*n2)/2
      s=sqrt((n1*n2*(n1+n2+1))/12)
      return((R-mu)/s)
    })
    lines(xi,z, col=j, ty="b")
    return(z)
  })  
})
names(WRSTidnn)=nms

WRSTidnnres=cbind(colMaxs(abs(WRSTidnn[[1]]), value = T),xi[colMaxs(abs(WRSTidnn[[1]]))],
                  colMaxs(abs(WRSTidnn[[2]]), value = T),xi[colMaxs(abs(WRSTidnn[[2]]))])

WRSTidnnBS=lapply(nms, function(iso){
  M=reabydiaBS[[iso]]
  X=pbsapply(M, function(m){
    sapply(1:4, function(j){
      z=sapply(xi, function(i){
        fc=cut(idnn, c(min(idnn)-1,i,max(idnn)+1),labels=F)
        Y1=m[fc==2,j]
        Y2=m[fc==1,j]
        res=wilcox.test(x=Y1,y=Y2)
        R=res$statistic
        n1=length(Y1)
        n2=length(Y2)
        mu=(n1*n2)/2
        s=sqrt((n1*n2*(n1+n2+1))/12)
        return((R-mu)/s)
      })
      return(c(xi[which.max(abs(z))], z[which.max(abs(z))]))
    })  
  })
  array(X, dim=c(2,4,100))
})
names(WRSTidnnBS)=nms

WidBS=lapply(WRSTidnnBS, function(Tz){
  t(array(Tz, dim=c(2,4*100)))
})

## Fig Idiotypy --------

spearCor=sapply(nms, function(iso) cor(reabydia[[iso]],idnn[rownames(reabydia[[iso]])], method = "spearman"))
rownames(spearCor)=colnames(reabydia$IgM)

plot(0:29,(WRSTidnn[[1]][1:30,1]), ty="l", ylim=c(-7,7),  col=2, yaxt="n",
     main="IgM Reactivity", xlab="Split point - # of NN idiotopes", xaxt="n",
     ylab="Wilcoxon RST Z", frame.plot=F)
axis(1,at=seq(0,30,2), labels=seq(0,30,2), pos=-7)
axis(2,at=seq(-7,7,2), labels=seq(-7,7,2), pos=0)
lines(0:29,(WRSTidnn[[1]][1:30,2]), pch=16, col=1)
lines(0:29,(WRSTidnn[[1]][1:30,3]), pch=16, col=3, lwd=2)
lines(0:29,(WRSTidnn[[1]][1:30,4]), pch=16,  col=4, lwd=2)
polygon(c(0,30,30,0), c(-1.96,-1.96,1.96,1.96), col=rgb(0.3,0.3,0.3,0.2))
legend(24,-4, legend=c("C","A","F","ND"), col=c(1,2,3,4), bty="n", lwd=2)
par(new=T, mai=c(6.2,3.5,1.2,1.2))
barplot(spearCor[c(2,1,3,4),1], col=1:4, ylim=c(-0.5,0.5), 
        ylab=expression(Spearman~rho), yaxt="n", xaxt="n")
axis(1, at=1:4, labels=c("C","A","F","ND"), cex.axis=0.8)
axis(2, at=c(-0.5,0,0.5), cex.axis=0.8)

plot(0:29,(WRSTidnn[[2]][1:30,1]), ty="l", ylim=c(-7,7),  col=2, yaxt="n",
     main="IgG Reactivity", xlab="Split point - # of NN idiotopes", xaxt="n",
     ylab="Wilcoxon RST Z", frame.plot=F)
axis(1,at=seq(0,30,2), labels=seq(0,30,2), pos=-7)
axis(2,at=seq(-7,7,2), labels=seq(-7,7,2), pos=0)
lines(0:29,(WRSTidnn[[2]][1:30,2]), pch=16, col=1)
lines(0:29,(WRSTidnn[[2]][1:30,3]), pch=16, col=3, lwd=2)
lines(0:29,(WRSTidnn[[2]][1:30,4]), pch=16,  col=4, lwd=2)
polygon(c(0,30,30,0), c(-1.96,-1.96,1.96,1.96), col=rgb(0.3,0.3,0.3,0.2))
legend(24,-4, legend=c("C","A","F","ND"), col=c(1,2,3,4), bty="n", lwd=2)
par(new=T, mai=c(6.2,3.5,1.2,1.2))
barplot(spearCor[c(2,1,3,4),2], col=1:4, ylim=c(-0.5,0.5), 
        ylab=expression(Spearman~rho), yaxt="n", xaxt="n")
axis(1, at=1:4, labels=c("C","A","F","ND"), cex.axis=0.8)
axis(2, at=c(-0.5,0,0.5), cex.axis=0.8)


heatmap.2()

# Clustering Leiden -------------------------------------------------------

clst=lapply(Gi, function(g){
  X=pbsapply(1:1000, function(j){
    cluster_leiden(g, objective_function = "modularity",resolution_parameter = 3)$membership
  })
  kmn=median(colMaxs(X, value = T))
  D=pbapply(X,1,function(l1){
    apply(X,1, function(l2){
      sum(l1==l2)
    })
  })
  D=log10(1/(D+0.5))
  D=(D-min(D))/diff(range(D))
  Xcl=hclust(as.dist(D), method="ward.D2")
  cutree(Xcl,k=kmn)
})
names(clst)=nms

clstmmb=clst
clst=lapply(nms, function(iso){
  L=clst[[iso]]
  x=aggregate(names(V(Gi[[iso]])), by=list(L), list)
  nm=x$Group.1;x=x$x; names(x)=nm
  return(x)
})
names(clst)=nms
clstN=sapply(clstmmb, function(cl) length(table(cl)))

#### Cluster bootstrap ------

clstBS=lapply(GiBS, function(giso){
  plan("multisession", workers=10)
  clmm=pblapply(giso, function(g){
    X=future_sapply(1:1000, function(j){
      X=cluster_leiden(g, objective_function = "modularity",resolution_parameter = 3)$membership
    }, future.seed=T)
    kmn=median(colMaxs(X, value = T))
    D=apply(X,1,function(l1){
      apply(X,1,function(l2){
        sum(l1==l2)
      })
    })
    D=log10(1/(D+0.5))
    D=(D-min(D))/diff(range(D))
    Xcl=hclust(as.dist(D), method="ward.D2")
    cutree(Xcl,k=kmn)
  })
  closeAllConnections()
  return(clmm)
})

clstBSns=lapply(clstBS, function(cliso){
  sapply(cliso, max)
})

clstBSmod=lapply(nms, function(iso){
  sapply(1:100, function(i){
    modularity(GiBS[[iso]][[i]], clstBS[[iso]][[i]])
  })
})

names(clstBSmod)=nms

clstMod=sapply(nms, function(iso) modularity(Gi[[iso]], clstmmb[[iso]]))

zclstNs=sapply(1:2, function(i){
  (clstN[i]-mean(clstBSns[[i]]))/sd(clstBSns[[i]])
})

zclstmods=sapply(1:2, function(i){
  (clstMod[i]-mean(clstBSmod[[i]]))/sd(clstBSmod[[i]])
})




#### Clustering IgM/IgG comparison ------------------------------------------------------------

clmx_=sapply(clst$IgM[lengths(clst$IgM)>1], function(c1){
  sapply(clst$IgG[lengths(clst$IgG)>1], function(c2){
    length(intersect(c1,c2))
  })
})
isng=rowSums(clmx_)>0
jsng=colSums(clmx_)>0
clmx_=clmx_[isng,jsng]

clmx=clmx_/sum(clmx_)
dimnames(clmx)=list(paste("C",1:nrow(clmx),"G", sep="_"),paste("C",1:ncol(clmx),"M", sep="_"))
X=lapply(clst, function(L){
  x=melt(L)
  xm=max(as.numeric(x$L1))
  ppmore=pp[!(pp %in% x$value)]
  addX=data.frame(ppmore,as.character(seq((xm+1), by=1, length.out=length(ppmore))))
  colnames(addX)=colnames(x)
  x=rbind(x,addX)
  return(x)
})

x=t(sapply(pp, function(p) {
  clim=sum(X$IgM$L1==X$IgM$L1[X$IgM$value==p])
  clig=sum(X$IgG$L1==X$IgG$L1[X$IgG$value==p])
  return(c(clim,clig))
}))

poff=rownames(x)[apply(x<10,1,all)]
Xno1=lapply(X, function(L){
  L[!(L$value %in% poff),]
})

amiMG=AMI(X$IgM$L1,X$IgG$L1)
amiMGno1=AMI(Xno1$IgM$L1,Xno1$IgG$L1)

heatmap.2(clmx_, hclustfun  = function(x) hclust(x, method = "ward.D2"), 
          dendrogram = "none", trace = "none", col=cpl(100), keysize = 1,
          density.info = "none", key.title = "", key.par = list(cex=0.5))

#### Cluster / diagnosis stats ---- 

clComparison=lapply(nms, function(iso){
  re=lapply(clst[[iso]], function(l){
    if(length(l)>1){
      n1=DnGM[[iso]][l,]
      n1=melt(n1)
      x=unlist(stri_extract_all(n1$Var2, regex="\\w+(?=_)"))
      n1$Var2=x
      tbl=dunn.test(n1$value, n1$Var2, method = "bh", kw=T)
      p=tbl$P
      res=tbl$Z
      nm=tbl$comparisons
      res=c(res,p)
      names(res)=c(nm,paste(nm,"p",sep="_"))
      return(res)
    } else return(NULL)
  })
  re=re[lengths(re)>0]
  re=t(as.data.frame(re))
  re[order(apply(abs(re[,1:6]), 1, mean), decreasing = T),]
})

clcomP=lapply(1:2, function(i) 
  array(p.adjust(clComparison[[i]][,7:12], method = "BH")<0.01, dim=dim(clComparison[[i]][,7:12])))
clComparison=lapply(1:2, function(i) {
  X=clComparison[[i]][rowSums(clcomP[[i]])>0,1:6]
  rownames(X)=paste(rownames(X), nms[[i]], sep="_")
  return(X)
})
names(clComparison)=nms

hclVcl=hclust(dist(rbind(clComparison$IgM,clComparison$IgG)), method="ward.D2")
plot(hclVcl)
clCompMG=rbind(clComparison$IgM,clComparison$IgG)[hclVcl$order,]
rownames(clCompMG)=sub("X","C", rownames(clCompMG))
isosupl=unlist(stri_extract_all(rownames(clCompMG), regex="(?<=_)\\w+"))

###### Fig biplot clusters -----

pcaclv=PCA(clCompMG)
fviz_pca_biplot(pcaclv, repel=T,palette=palette()[c(4,7)], col.ind=as.factor(isosupl), axes = c(1,2), col.var="red")
fviz_pca_biplot(pcaclv, repel=T,palette=palette()[c(4,7)], col.ind=as.factor(isosupl), axes = c(1,3), col.var="red")
fviz_screeplot(pcaclv)

#### Graph bootstrap characteristics report ------------------------------------------
###### Fig 4charact  -----------------------
par0=par()
par(mfrow=c(1,4))
par(bty="n")
par(mai=par("mai")+c(0,0.1,0,0))

yr=c(0.015, 0.04)
boxplot(bsdns, notch=T, ylim=yr,xlim=c(0,2.5), ylab="Graph density", 
        yaxt="n", xaxt="n", cex.lab=1.5)
axis(1,at=c(0,1,2,2.5),labels=c("","IgM","IgG",""), pos=0.015, cex.axis=1.5)
axis(2,at=seq(from=0.015, to=0.04, by=0.005),pos=0, cex.axis=1.2)
points(c(graph.density(Gi$IgM),graph.density(Gi$IgG)), col=2, pch=16, cex=1.5)
text(c(1.15,1.95), c(0.039,0.033), labels=paste("z=",round(zdns,2)), cex=1.2)


yr=range(unlist(c(compBSlcm,grcomp)))+c(-8,11)
boxplot(compBSlcm,notch=T, ylim=yr,xlim=c(0,2.5),
        ylab="Size of large component", yaxt="n", xaxt="n", cex.lab=1.5)
axis(1,at=c(0,1,2,2.5),labels = c("","IgM","IgG",""),pos=660, cex.axis=1.5)
axis(2,at=seq(from=yr[1], to=yr[2], by=10), pos=0, cex.axis=1.2)
points(grcomp, col=2, pch=16, cex=1.5)
text(c(1.15,2.02), c(711, 703), labels=paste("z=",round(zcomp,2)), cex=1.2)


yr=c(0,100)
boxplot(clstBSns,notch=T, ylim=yr, xlim=c(0,2.5), 
        ylab="Number of clusters", yaxt="n", xaxt="n", cex.lab=1.5)
axis(1,at=c(0,1,2,2.5),labels = c("","IgM","IgG",""), pos=yr[1], cex.axis=1.5)
axis(2,at=seq(from=yr[1], to=yr[2], by=floor(diff(yr)/5)), pos=0, cex.axis=1.2)
points(clstN, col=2, pch=16, cex=1.5)
text(c(1.25,2), clstN-2, labels=paste("z=",round(zclstNs,2)), cex=1.2)


yr=c(0.15,0.30)
boxplot(clstBSmod,notch=T, ylim=yr, xlim=c(0,2.5), cex.lab=1.5, 
        ylab="Modularity of clustering", xaxt="n", yaxt="n")
axis(1,at=c(0,1,2,2.5),labels = c("","IgM","IgG",""), pos=yr[1], cex.axis=1.5)
axis(2,at=seq(from=yr[1], to=yr[2], by=floor(100*diff(yr)/5)/100), pos=0, cex.axis=1.2)
points(clstMod, col=2, pch=16, cex=1.5)
text(c(1.15,2), clstMod+c(0.002,-0.002), labels=paste("z=",round(zclstmods,2)), cex=1.2)


###### Fig cliques ---------------------


X=lapply(maxclqBS, function(L) melt(L))

par(mfrow=c(1,2))
par(bty="n")
for (iso in nms){
  n=barplot(table(lengths(maxclq[[iso]])), main=iso, xlab="Clique size", 
            ylab="Counts", xlim=c(1,9), ylim=c(0,3500), col="lightgrey")
  boxplot(X[[iso]]$value~X[[iso]]$Var1,  xlim=c(1,9), ylim=c(220,6000), add=T,
          at=n[seq_along(unique(X[[iso]]$Var1))],
          col=rgb(0,0,0,0.5), xaxt="n", yaxt="n", xlab="", ylab="", cex=0.3)
  lines(aggregate(X[[iso]]$value, by=list(X[[iso]]$Var1), median))
  par(new=F)
}

###### Fig Degree ---------

deGi=sapply(Gi, degree)

deGiBS=lapply(GiBS, function(L){
  sapply(L,degree)
})

h=hist(log10(unlist(deGiBS)+0.5))
br=c(-0.5,0,h$breaks[-(1:9)])
h=hist(log10(unlist(deGiBS)+0.5), breaks=br)
deGiBSt=lapply(deGiBS,function(L){
  apply(L,2,function(l){
    h=hist(log10(l+0.5), breaks=br, plot=F)
    return(h$counts)
  })
})

hmids=h$mids
deGit=sapply(deGi,function(cl) hist(log10(cl+0.5), breaks=br, plot=F)$counts)

par(mfrow=c(1,2))
par(bty="n")
for (iso in nms){
  X=deGit[,iso]+0.5
  hi=barplot(X, main=iso, xlab="Degree", ylab="Counts", pch=16,  
           xaxt="n", yaxt="n", ylim=c(0.5,115.5), col=rgb(1,0,0,0.3))
  axis(1, at=hi, labels=floor(10^(h$mids)), las=2)
  axis(2, at=seq(0.5, 115.5,5), labels=seq(0, 115,5))#
  boxplot(t(deGiBSt[[iso]]),at=hi,  ylim=c(0.5,115.5), add=T,
          col=rgb(0,0,0,0.3), xaxt="n", yaxt="n", xlab="", ylab="", cex=0.3)
  lines(hi,apply(deGiBSt[[iso]], 1, median))
  par(new=F)
}



# Edge clustering --------

EE=lapply(nms, function(iso){
  ends(Gi[[iso]], E(Gi[[iso]]))
})
names(EE)=nms
EEclust=lapply(nms, function(iso){
  Eiso=EE[[iso]]
  D=t(apply(Eiso, 1, function(l){
   colMeans(DnGM[[iso]][l,])
  }))
  D=Dist(D)
  diag(D)=0
  D[D>median(D[lower.tri(D)])]=0
  GEE=graph_from_adjacency_matrix(D, mode="lower", weighted=T)
  w=max(edge_attr(GEE)$weight)
  cutsc=pbsapply(seq(0.4,0.55,length.out=20), function(wi){
    g=delete.edges(GEE, E(GEE)[edge_attr(GEE)$weight>wi])
    cl=cluster_leiden(g, objective_function = "modularity")
    tcl=table(cl$membership)
    c(wi,log10(prod(tcl[tcl>5])))
  })
  ct=cutsc[1,which.max(cutsc[2,])]
  GEE=delete_edges(GEE, E(GEE)[edge_attr(GEE)$weight>ct])
  cluster_leiden(GEE, objective_function = "modularity")$membership
})
names(EEclust)=nms
Gi=lapply(c("IgM","IgG"), function(iso) {
  G=Gi[[iso]]
  cm=components(G)
  G=induced.subgraph(G, V(G)[cm$membership==as.numeric(names(table(cm$membership)))[which.max(table(cm$membership))]])
  n=vcount(G)
  L=embed_laplacian_matrix(G,no=n-1, type = "I-DAD") 
  if (iso=="IgM") jj=705:708 else jj=696:699
  uMx=umap(L$X[,jj], min_dist = 0.005, spread=10,  n_epochs = 1500, 
           verbose=F, n_components = 2, init = "normlaplacian",n_sgd_threads=0)
  G=set_vertex_attr(G, name="x", value=70*uMx[,1])
  G=set_vertex_attr(G, name="y", value=70*uMx[,2])
  G=set_edge_attr(G, name="Ecl", value=EEclust[[iso]])
  write.graph(G, format = "graphml", file=paste("G",iso,".graphml", collapse="",sep=""))
  return(G)
}) 
names(Gi)=c("IgM","IgG")

# Combined clustering ------
cpl3=colorRampPalette(c("#000000FF","#00FFFFFF","#10FA10FF","#FFFF009F","#FFA0009F","#B50000"), alpha=T)
cpl5=c("#000000FF","#00FFFFFF","#10FA10FF","#FFFF009F","#FFA0009F","#B50000")

for (iso in nms) colnames(EE[[iso]])=c("V1","V2")

vEEclust=lapply(nms, function(iso){
  CL=clst[[iso]]
  lngee=unique(melt(data.frame(
    EE[[iso]],Cl=edge_attr(Gi[[iso]])$Ecl, stringsAsFactors = F),id.vars="Cl")[,-2])
  X=lapply(CL, function(cl){
    if (length(cl)>2) {
      g=induced.subgraph(Gi[[iso]], V(Gi[[iso]])[cl])
      ecl=edge_attr(g)$Ecl
      L=lngee[lngee[,1] %in% ecl & lngee[,2] %in% cl,1:2]
      x=aggregate(L$value, by=list(L$Cl), c)
      nm=x$Group.1; x=x$x; names(x)=nm
      return(x)
    } else return(NULL)
  })
  X=unlist(X, recursive=F)
  X[lengths(X)>2]
})
names(vEEclust)=c("IgM","IgG")  

clijstat=lapply(1:2, function(j){
  clij=vEEclust[[j]]
  if (j==1) M=DnGM$IgM else M=DnGM$IgG
  X=t(sapply(clij, function(L){
    mL=melt(M[L,])
    gr=rep(diapat$dia, eac=length(L))
    dun=dunn.test(mL$value,gr, method="none")
    Z=dun$Z; names(Z)=dun$comparisons
    P=dun$P; names(P)=paste(names(Z),"_p",sep="_")
    return(c(Z,P))
  }))
  rownames(X)=names(clij)
  X[,7:12]=p.adjust(X[,7:12], method="BH")
  X[,1:6]=X[,1:6]*(X[,7:12]<1e-3)
  X[,7:12]=X[,7:12]*(X[,7:12]<1e-3)
  X=X[rowSums(X)!=0,]
  clX=hclust(dist(X[,1:6]))
  x=cutree(clX, h=7)
  X=X[clX$order,]
  return(X)
})
names(clijstat)=nms

X=clijstat$IgM[,1:6]
X[,3]=-X[,3]
colnames(X)[3]="F-C"
X[,5]=-X[,5]
colnames(X)[5]="ND-C"
X=PCA(X)
fviz_screeplot(X)
fviz_pca_biplot(X,  repel=T, col.var="red", col.ind = "grey")

X=clijstat$IgG[,1:6]
X[,3]=-X[,3]
colnames(X)[3]="F-C"
X[,5]=-X[,5]
colnames(X)[5]="ND-C"
X=PCA(X)
fviz_screeplot(X)
fviz_pca_biplot(X, repel=T, col.var="red", col.ind = "grey")
fviz_pca_biplot(X, repel=T, col.var="red", col.ind = "grey", axes=c(1,3))

bestMclij=clijstat$IgM [rowsums(abs(clijstat$IgM)>5)>0,]
hclij=hclust(dist(bestMclij), method = "ward.D2")
plot(hclij)
clfinalM=cutree(hclij, h=5)
ppclfinalM=melt(vEEclust$IgM[names(clfinalM)])
ppclfinalM=cbind(ppclfinalM,cl=clfinalM[ppclfinalM$L1])
ppclfinalM=aggregate(ppclfinalM$value, by=list(ppclfinalM$cl), function(x) unique(unlist(x)))
nm=c("C>ND.F.A","C>F.A&ND>A","C>F.A","ND.F>C","A.F>ND.C","A>C","A.F>C")
ppclfinalM=ppclfinalM$x
names(ppclfinalM)=nm

bestGclij=clijstat$IgG[rowsums(abs(clijstat$IgG)>5)>0,]
hclij=hclust(dist(bestGclij), method = "ward.D2")
plot(hclij)
clfinalG=cutree(hclij, h=3)
ppclfinalG=melt(vEEclust$IgG[names(clfinalG)])
ppclfinalG=cbind(ppclfinalG,cl=clfinalG[ppclfinalG$L1])
ppclfinalG=aggregate(ppclfinalG$value, by=list(ppclfinalG$cl), function(x) unique(unlist(x)))
nm=c("A.C>F","A>ND","F.ND>C","F>C")
ppclfinalG=ppclfinalG$x
names(ppclfinalG)=nm

ppclfinalM=lapply(ppclfinalM, unique)
write.table(sapply(ppclfinalM, function(l) c(l,rep(" ",max(lengths(ppclfinalM))-length(l)))), file="ppclfinalM.txt")
ppclfinalG=lapply(ppclfinalG, unique)
write.table(sapply(ppclfinalG, function(l) c(l,rep(" ",max(lengths(ppclfinalG))-length(l)))), file="ppclfinalG.txt")

ppclfinMId=sapply(ppclfinalM, function(L) {
  log10(idnn[L]+1)
})
ppclfinGId=sapply(ppclfinalG, function(L) {
  log10(idnn[L]+1)
})

ppclfinalMId=c(ppclfinMId,All=list(log10(idnn[!(idnn %in% unlist(ppclfinalM))]+1)))
ppclfinalMId=ppclfinalMId[order(sapply(ppclfinalMId, median))]
par(mai=c(1.5,0.82,0.82,0.42))
vioplot(ppclfinalMId, h=0.15, ylab="N idiotope neighbors", ylim=c(0,2.5), frame.plot=F, xaxt="n", yaxt="n")
axis(1,at=seq(from=0, to=9, by=1), labels = c("",names(ppclfinalMId),""), pos=0,las=2)
axis(2,at=seq(from=0, to=2.5, by=0.5))
text(c(3.5, 5.5),c(1.5, 1.5),labels = c("*","*"), cex=2, col=c("red","black"))
dt=dunn.test(ppclfinalMId)

ppclfinalGId=c(ppclfinGId,All=list(log10(idnn[!(idnn %in% unlist(ppclfinalG))]+1)))
ppclfinalGId=ppclfinalGId[order(sapply(ppclfinalGId, median))]
par(mai=c(1.5,0.82,0.82,0.42))
vioplot(ppclfinalGId, h=0.15, ylab="N idiotope neighbors", ylim=c(0,2.5), frame.plot=F, xaxt="n", yaxt="n")
axis(1,at=seq(from=0, to=6, by=1), labels = c("",names(ppclfinalGId),""), pos=0,las=2)
axis(2,at=seq(from=0, to=2.5, by=0.5))
lines(c(3,4),c(2.45,2.45));lines(c(3,3),c(2.4,2.45));lines(c(4,4),c(2.4,2.45))
text(3.5,2.5, labels="*")
dt=dunn.test(ppclfinalGId)

pM=unique(unlist(ppclfinalM))
pG=unique(unlist(ppclfinalG))
length(pM)
length(pG)
length(intersect(pM,pG))

ppfinal=c(ppclfinalM,ppclfinalG)
names(ppfinal)=paste(names(ppfinal), c(rep("IgM",length(ppclfinalM)),rep("IgG",length(ppclfinalG))),sep=":")
ppfinal_all=c(ppfinal,AllIgM=list(pp[!(pp %in% unique(unlist(ppfinal)))]),
              AllIgG=list(pp[!(pp %in% unique(unlist(ppfinal)))]))

corbydiaj=sapply(seq_along(ppfinal_all), function(i){
  x=ppfinal_all[[i]]
  m=t(sapply(x, function(p){
    if (i %in% c(1:7,12)) X=DnGM$IgM[p,] else X=DnGM$IgG[p,]
    qs=quantile(X, c(0.25,0.5,0.75))
    m=(X-qs[2])/(diff(qs[c(1,3)]))
  }))
  for (p in x){
    plot(m[p,order(diapat$dia)], ty="p", ylim=range(m), 
         col=rep(c("#FF000060","#00000060","#0000FF60","#00FF0060"), each=5), 
         main=names(ppfinal_all)[i], pch=16)
    lines(1:20,m[p,order(diapat$dia)], col=rgb(0,0,0,0.3))
    par(new=T)
  }
  par(new=F)
  sapply(unique(diapat$dia), function(d){
    cm=(1/dist(t(m[,diapat$dia==d]), method = "manhattan"))
    return(mean(cm))
  })
})

colnames(corbydiaj)=names(ppfinal_all)

mnbydiaj=sapply(seq_along(ppfinal_all), function(i){
  x=ppfinal_all[[i]]
  if (i %in% c(1:7,12)) X=DnGM$IgM[x,] else X=DnGM$IgG[x,]
  sapply(unique(diapat$dia), function(d){
    return(mean(X[,diapat$dia==d]))
  })
})

colnames(mnbydiaj)=names(ppfinal_all)

corrplot(cor(mnbydiaj,corbydiaj), method = "color", col=cpl(100), ylab="Reactivtiy")
      #,order = "hclust", hclust.method = "ward.D2")
corrplot(cor(corbydiaj), method = "color", col=cpl(100), order = "hclust")
corrplot(cor(mnbydiaj), method = "color", col=cpl(100), is.corr = F, order = "hclust")

## Fig combined clustering -----

X=clst$IgM[[6]]

gx=induced_subgraph(Gi$IgM, V(Gi$IgM)[X])
n=vcount(gx)
L=embed_laplacian_matrix(gx,no=n-1, type = "I-DAD")
plot(L$D)
uL=umap(L$X[,(n-3):(n-1)],  n_neighbors = 4, min_dist = 0.005, spread=10,  n_epochs = 1500,
    verbose=F, n_components = 2, init = "normlaplacian",n_sgd_threads=1, seed=878234)
V(gx)$label.cex=1
ecol=cut(edge_attr(gx)$Ecl, length(unique(edge_attr(gx)$Ecl)), labels=F)
plot(gx, layout=uL, vertex.size=0, edge.color=ecol, edge.width=3)
gx=set_edge_attr(gx, name="ecol", value=ecol)
write_graph(gx, format="graphml", file="gx.graphml")


ppCMlow=unique(unlist(ppfinal[4:7]))
ppAFG=unique(unlist(ppfinal[8:11]))
isosw=intersect(ppCMlow,ppAFG)
chisq.posthoc.test(table(names(idnn) %in% isosw,idnn>2), simulate.p.value=T)


# BLAST ----

getProteome(organism="9606", release=76)

# Drost HG, Paszkowski J. Biomartr: genomic data retrieval with R. Bioinformatics (2017) 33(8): 1216-1217. doi:10.1093/bioinformatics/btw821.
# "_ncbi_downloads/proteomes/9606_protein_refseq.faa.gz" renamed to HuProt.faa

# BLAST arguments are form https://www.ncbi.nlm.nih.gov/books/NBK279690/

HuProt=seqinr::read.fasta(file="HuProt.faa")
HuProtannot=lapply(HuProt, function(x) attributes(x)$Annot)
zfngr=grep("zinc finger", HuProtannot)
length(zfngr)/length(HuProtannot)

bl=blast(db="HuProt.faa",remote=F, type = "blastp")

BLall=pblapply(pp, function(ppL) {
  X=predict(bl, AAStringSet(ppL), 
            BLAST_args=list("-word_size 2", "-gapopen 9","-gapextend 1", 
                            "-matrix PAM30", "-threshold 16", 
                            "-comp_based_stats 0", "-window_size 15"))
  X[X$mismatch<2,]
})
names(BLall)=pp

j=sapply(BLall,nrow); sum(j>0)/length(pp)

BLALL=BLall[j>0]

allprotint=unique(unlist(sapply(BLALL,function(x) x$sseqid)))
allprotint=HuProt[allprotint]
allprotintann=unique(lapply(allprotint, function(x) attributes(x)$Annot))
zfngrint=grep("zinc finger", allprotintann)
length(zfngrint)/length(allprotint)
grep("amyloid",allprotintann, value=T)
grep("protein tau",allprotintann, value=T)
grep("glial fibrillary",allprotintann, value=T)
grep("chitinase-3-like protein 1",allprotintann, value=T)
grep("triggering receptor expressed on myeloid cells",allprotintann, value=T)
grep("galectin-3",allprotintann, value=T)
grep("alpha-2-macroglobulin",allprotintann, value=T)  # tolerance?


distrblall=sapply(BLALL, nrow)
hist(log10(distrblall))
BLALL[distrblall>100]
plot(distrblall,idnn[names(distrblall)], log="xy")
summary(lm(distrblall~idnn[names(distrblall)]))
plot(aggregate(distrblall, by=list(idnn[names(distrblall)]), "mean"), ty="l")
wilcox.test(distrblall[ppf],distrblall[!(names(distrblall) %in% ppf)])

ppf=unique(unlist(ppfinal))
ppfbl=ppf[ppf %in% names(BLALL)]
pptHits=lapply(ppfbl, function(p) { 
    L=BLALL[[p]]
    X=t(sapply(seq_along(L$sseqid), function(i) {
        X=L$sseqid[i]
        x=attributes(HuProt[[X]])$Annot
        x2=gsub("[ ]$","",gsub("\\[Homo sapiens]|>[^ ]*|\\ isoform \\w*|\\ member \\w*| precursor \\w*","",x))   
        S=toupper(paste(HuProt[[X]][L$sstart[i]:L$send[i]], collapse=""))        
        c(X,x2,S,p)
      }))
      X=aggregate(list(X[,1]), by=list(X[,2],X[,3],X[,4]), "list")
      colnames(X)=c("Protein","Protein Seq","Mimotope Seq", "IDs")
      j=nchar(X$`Protein Seq`)>=6
      return(X[j,])
})
names(pptHits)=ppfbl

x=pptHits[[1]]
for (i in 2:length(pptHits)){
  if (length(pptHits[[i]])>0) x=rbind(x,pptHits[[i]])
}

pptHitsall=x
ppthitsbycl=lapply(ppfinal,function(l) {
  x=as.data.frame(pptHitsall[pptHitsall[,"Mimotope Seq"] %in% l,], )
  x$IDs=sapply(x$IDs, function(l) paste(l,collapse="|"))
  return(x)
})


ppthitsbycl=melt(ppthitsbycl)
write.table(ppthitsbycl,file="ppHits_bycl.txt")



tx=rbind(c(length(BLALL),length(pp)-length(BLALL)),
         c(length(ppfbl), length(ppfinal_u)-length(ppfbl)))

chisq.test(tx)

chisq.test(table(idnn>0,pp %in% ppfinal_u))

zfngrHits=length(grep("zinc finger",pptHits))
zfngrHits/length(pptHits)


ppfinal_tu=sapply(ppf, function(cl) ppf %in% cl)
colnames(ppfinal_tu)=sub("\\&"," ",colnames(ppfinal_tu))
rownames(ppfinal_tu)=ppfinal_u
vnn=euler(ppfinal_tu,shape="ellipse", input="disjoint")
plot(vnn, quantities=T)

mean(aggregate(pptHits$Group.1, by=list(pptHits$L1), "length")[,2])

ppthitsbyclbyprot=aggregate(ppthitsbycl$L1, by=list(ppthitsbycl$Protein), "list")


justpphits=unique(ppthitsbycl$`Mimotope Seq`)


MM=DnGM$IgM[justpphits,]; MG=DnGM$IgG[justpphits,]
Mx=t(cbind(scale(t(MM)),scale(t(MG))))

plot3d(cmdscale(dist(t(DnGM$IgG[dd,])),k=3), col=diall, size=10)
rownames(Mx)=paste(rownames(Mx), rep(c("M","G"),each=length(justpphits)), sep="_")
dd=rownames(rfe(Mx, as.factor(diall))[[1]])
dd=rownames(rfe(Mx[,diall %in% c(1,3)], as.factor(diall[diall %in% c(1,3)]))[[1]])
ddneat=cbind(unlist(stri_extract_all(dd, regex="\\w+(?=_)")),dd)
plot3d(cmdscale(dist(t(Mx[dd,])),k=3), col=diall, size=10, xlab="D1", ylab="D2", zlab = "D3")

plot3d(cmdscale(dist(t(Mx[dd,])),k=3), col=diall, size=10, xlab="D1", ylab="D2", zlab = "D3")
plot(cmdscale(dist(t(Mx[dd,])),k=2), pch=16, col=diall, cex=2, xlab="D1", ylab="D2")

write.csv(ddneat, file="ddneat.csv")

