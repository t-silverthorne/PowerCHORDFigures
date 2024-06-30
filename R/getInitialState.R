getInitialState=function(freq,Nfine,Nmeas){
  Cmats = getFourQuadBlocks(freq,Nfine,Nmeas)
  mu = rep(0,Nfine)
  mu[sample(1:Nfine,Nmeas)]=1

  a_0  = t(mu)%*%Cmats$Cm11%*%mu 
  b_0  = t(mu)%*%Cmats$Cm12%*%mu 
  c_0  = t(mu)%*%Cmats$Cm22%*%mu 
  p1_0 = sqrt(a_0^2 - 2*a_0*c_0 + 4*b_0^2 + c_0^2)
  p2_0 = a_0/2 - c_0/2 - p1_0/2
  p3_0 = a_0*p2_0 
  p4_0 = b_0^2 
  p5_0 = p3_0*p2_0 +p4_0*(c_0+2*p2_0) 
  p6_0 = p2_0^2 + b_0^2
  p7_0 = 1/p6_0
  return(c(mu,a_0,b_0,c_0,p1_0,p2_0,p3_0,p4_0,p5_0,p6_0,p7_0))
}