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
  