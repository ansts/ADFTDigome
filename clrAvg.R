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