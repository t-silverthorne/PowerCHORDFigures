getPcQuadCsts=function(xfreq,fmin,fmax,Nfreq,Nfine,Nmeas){
  QC=list() # list of quadratic constraints
  
  # matrices for regularization term 
  tau  = c(1:Nfine)/Nfine -1/Nfine 
  fvec = seq(from=fmin,to=fmax,length.out=Nfreq)
  
  for (ii in c(1:Nfreq)){
    xf     = xfreq[ii] 
    freq   = fvec[ii]
    cvec   = cos(2*pi*freq*tau)
    svec   = sin(2*pi*freq*tau)
    
    C1top  = cbind(diag(cvec*cvec),diag(cvec*svec))
    C1bot  = cbind(diag(cvec*svec),diag(svec*svec))
    C1     = rbind(C1top,C1bot)
     
    C2top  = cbind(outer(cvec,cvec,'*'),outer(cvec,svec,'*'))
    C2bot  = cbind(outer(cvec,svec,'*'),outer(svec,svec,'*'))
    C2     = rbind(C2top,C2bot)
     
    Cmat   = C1-C2/Nmeas
     
    Cm11   = Cmat[1:Nfine,1:Nfine]
    Cm12   = Cmat[1:Nfine,(Nfine+1):(2*Nfine)]
    Cm21   = Cmat[(Nfine+1):(2*Nfine),1:Nfine]
    Cm22   = Cmat[(Nfine+1):(2*Nfine),(Nfine+1):(2*Nfine)]
   
    qcmat  = xf*Cm11 + (1-xf)*Cm22 + sqrt(xf)*sqrt(1-xf)*(Cm12+Cm21)
    qcmat  = cbind(rbind(qcmat,rep(0,Nfine)),rep(0,Nfine+1))
    qc     = c(rep(0,Nfine),-1)
    beta   = 0
    sense  = '>'
    QC[[ii]] = list(Qc=qcmat,q=qc,rhs=beta,sense=sense)
  }
  return(QC)
} 