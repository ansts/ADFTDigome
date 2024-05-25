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