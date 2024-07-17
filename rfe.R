
# A function which iteratively scans the set of features and finds at each cycle the feature  
# the removal of which imrpoves most the separation between the classes of interest c, then removes it 
# leaving one less feature for the next cycle. The separation is measured using a composite clustering criterion 
# consisting of the sum of Dunn's, BHgamma criteria as well as negative Connectivity.
# The criteria are first transformed to compensate for their dependence on the dimensionality 
# and to restrict them in the range between -1 and 1. The conversion functions dunnfix, connfix and bhgamfix 
# are vectorized for faster calculations. 

rfe=function(m0,c){
  require(clValid)
  require(parallel)
  require(cluster)
  require(Rfast)
  require(reshape2)
  require(matrixStats)

  if (is.logical(c)) cc=as.integer(c)+1 else cc=as.integer(c)
  #pt=proc.time()
  D0=c()
  Ds0=c()
  bad=list()
  bads=list()
  m=m0
  N0=nrow(m0)
  nm=rownames(m0)
  nc=unique(c)
  vconf=Vectorize(connfix)
  vbhgf=Vectorize(bhgamfix)
  vdunf=Vectorize(dunnfix)
  repeat {
      D=c()
      N=nrow(m)
      if (N<3) break
      D=t(sapply(1:N, function(i) {
            d1=Dist(t(m[-i,]))
            x=c(BHgamma(d1,cc), dunn(as.dist(d1),cc),connectivity(as.dist(d1),cc))
            return(x)
        }))
      NN=rep(N-1,N)
      D=cbind(vbhgf(D[,1]),vdunf(D[,2],NN),-vconf(D[,3]))
      D=rowsums(D)
      D=c(which.max(D), max(D))                       
      
      md=D[[1]]                    
      dmin=D[[2]]
      bad=c(bad,unlist(rownames(m)[md]))
      m=m[-md,]
      D0=rbind(D0,c(N-1,dmin)) 
      # print(c(N,D[[2]]))  
      fnm="scanlog.txt"
      cat(c(N,D[[2]]), file=fnm, sep="\n", append=T)
  }
  
  bad=c(unlist(bad),unlist(rownames(m)))
  # yr=range(D0[,2])
  # plot(nrow(D0):1,D0[,2])


  if (length(D0)>0) crit=D0[which.max(D0[,2]),1] else return(list(NA,NA,NA))
  if (crit<2) crit=2
  m=m0[!(nm %in% bad[1:(length(bad)-crit)]),]
  
  return(list(m,bad,D0))
}

# Functions correcting the clustering criteria dependence on the dimensionality
# Valid for dimensions<4200

dunnfix=function(x,i){
  a_1=4.777408e-01
  b_1=1.312427e+01
  a_2=1.336025e-06 
  a0=5.507259e-02
  a1=-8.507095e-05
  a2=7.159110e-08
  a3=-3.028945e-11
  a4= 6.110248e-15
  a5= -4.700154e-19
  
  y=x-(i*a_1/(i+b_1)+i*a_2)
  y=y/(a5*i^5+a4*i^4+a3*i^3+a2*i^2+a1*i+a0)           
  return(y)
}

bhgamfix=function(x){
  y=(x-4.541866e-05)/0.08585348                                        #y=x*(i+342.277)/87.12
  return(y)
}
connfix=function(x){
  y=(x-30.74765)/3.571313                         
  return(y)                                  
}

# BH gamma clustering criterion implementation

BHgamma=function(A,c){
  n=ncol(A)
  A=as.dist(A)
  N=length(A)
  cc=matrix(rep(c,n),n)
  B=as.integer(as.dist(cc!=t(cc)))
  A=matrix(rep(A,N),N)
  sm=A-t(A)
  r=range(sm)
  eps=diff(r)/100
  sp=sm>eps
  sm=sm<(-eps)
  sp=sum(sp[B==1,B==0])
  sm=sum(sm[B==1,B==0])
  return((sp-sm)/(sp+sm))
}

