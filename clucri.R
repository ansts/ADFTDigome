
# A function which measures the cluster separation using a composite clustering criterion 
# consisting of the best value selected among Dunn's or BHgamma criteria as well as negative Connectivity.
# The criteria are transformed to compensate for their dependence on the dimensionality 
# and to restrict them in the range between -1 and 1. The conversion functions dunnfix, connfix and bhgamfix 
# are vectorized for faster calculations.

clucri=function(m0,c){
  require(clValid)
  require(cluster)
  require(Rfast)

  if (is.logical(c)) cc=as.integer(c)+1 else cc=as.integer(c)
 
  m=m0
  N0=nrow(m0)
  nm=rownames(m0)
  nc=unique(c)
  vconf=Vectorize(connfix)
  vbhgf=Vectorize(bhgamfix)
  vdunf=Vectorize(dunnfix)
  N=nrow(m)
  d1=Dist(t(m))
  x=c(BHgamma(d1,cc), dunn(as.dist(d1),cc),connectivity(d1,cc,neighbSize = min(length(cc) %/% 2,10)))
  D=sum(vbhgf(x[1]),vdunf(x[2],N),-vconf(x[3]))

  return(D)
}

bhgamfix=function(x){
  y=(x-4.541866e-05)/0.08585348                                        
  return(y)
}

connfix=function(x){
  y=(x-30.74765)/3.571313                         
  return(y)                                  
}

#Valid for dimensions<4200

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
