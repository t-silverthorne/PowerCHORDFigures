getFourQuadBlocks=function(freq,Nfine,Nmeas){
  # fine temporal grid
  tau    = c(1:Nfine)/Nfine -1/Nfine 
  cvec   = cos(2*pi*freq*tau)
  svec   = sin(2*pi*freq*tau)
  
  # build the main matrix
  C1top  = cbind(diag(cvec*cvec),diag(cvec*svec))
  C1bot  = cbind(diag(cvec*svec),diag(svec*svec))
  C1     = rbind(C1top,C1bot)
    
  C2top  = cbind(outer(cvec,cvec,'*'),outer(cvec,svec,'*'))
  C2bot  = cbind(outer(cvec,svec,'*'),outer(svec,svec,'*'))
  C2     = rbind(C2top,C2bot)
    
  Cmat   = C1-C2/Nmeas
    
  # split into four blocks of interest 
  Cm11   = Cmat[1:Nfine,1:Nfine]
  Cm12   = Cmat[1:Nfine,(Nfine+1):(2*Nfine)]
  Cm21   = Cmat[(Nfine+1):(2*Nfine),1:Nfine]
  Cm22   = Cmat[(Nfine+1):(2*Nfine),(Nfine+1):(2*Nfine)]

  return(list(Cm11=Cm11,Cm12=Cm12,Cm21=Cm21,Cm22=Cm22))
}