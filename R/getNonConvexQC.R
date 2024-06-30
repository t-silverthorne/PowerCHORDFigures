require(Matrix)
getNonConvexQC = function(Nfine,A11,A12,A22){
  # unpack
  n   = Nfine
  A11 = Matrix::.dense2sparse(A11)
  A12 = Matrix::.dense2sparse(A12)
  A22 = Matrix::.dense2sparse(A22)

  # indices for constraint forms
  ia   = 1 # matrix element top
  ib   = 2 # matrix element off-diag
  ic   = 3 # matrix element bottom
  ip1  = 4 # evec param 
  ip2  = 5 # evec param 
  ip3  = 6 # evec param 
  ip4  = 7 # evec param 
  ip5  = 8 # evec param 
  ip6  = 9 # evec param 
  ip7  = 10 # evec param 

  # number of pars
  np   = 10 
  QC   = list()

  # construct the discrete quadratic forms
  A11sp = bdiag(A11,Matrix(rep(0,np*np),nrow=np))
  A12sp = bdiag(A12,Matrix(rep(0,np*np),nrow=np))
  A22sp = bdiag(A22,Matrix(rep(0,np*np),nrow=np))

  q11   = sparseMatrix(1,j=ia,x=-1,dims=c(1,np))
  q12   = sparseMatrix(1,j=ib,x=-1,dims=c(1,np))
  q22   = sparseMatrix(1,j=ic,x=-1,dims=c(1,np))

  q11sp = cbind(Matrix(rep(0,n),nrow=1),q11) 
  q12sp = cbind(Matrix(rep(0,n),nrow=1),q12) 
  q22sp = cbind(Matrix(rep(0,n),nrow=1),q22) 

  QC[[1]]=list(Qc=A11sp,q=toSparseVector(q11sp),rhs=0,sense='=')
  QC[[2]]=list(Qc=A12sp,q=toSparseVector(q12sp),rhs=0,sense='=')
  QC[[3]]=list(Qc=A22sp,q=toSparseVector(q22sp),rhs=0,sense='=')

  # sqrt eigenvector term 
  Qp1     = sparseMatrix(ia,ia,x=1,dims=c(np,np))+
            sparseMatrix(ia,ic,x=-2,dims=c(np,np))+
            sparseMatrix(ib,ib,x=4,dims=c(np,np))+
            sparseMatrix(ic,ic,x=1,dims=c(np,np))+
            sparseMatrix(ip1,ip1,x=-1,dims=c(np,np))
  Qp1sp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qp1)
  QC[[4]] = list(Qc=Qp1sp,q=0,rhs=0,sense='=')

  #TODO: remove since now redundant
  #qp2     = sparseMatrix(c(1,1,1,1),j=c(ia,ic,ip1,ip2),
  #                        x=c(1/2,-1/2,-1/2,-1),dims=c(1,np))
  #qp2sp   = cbind(Matrix(rep(0,n),nrow=1),qp2) 
  #QC[[5]] = list(Qc=Matrix(rep(0,(Nfine+npar)^2),nrow=Nfine+npar),
  #          q=toSparseVector(qp2sp),rhs=0,sense='=')


  # 
  Qp3     = sparseMatrix(ia,ip2,x=1,dims=c(np,np))
  qp3     = sparseMatrix(c(1),j=c(ip3),x=c(-1),dims=c(1,np))
  Qp3sp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qp3)
  qp3sp   = cbind(Matrix(rep(0,n),nrow=1),qp3) 
  QC[[5]] = list(Qc=Qp3sp,q=toSparseVector(qp3sp),rhs=0,sense='=')


  # 
  Qp4     = sparseMatrix(ib,ib,x=1,dims=c(np,np))
  qp4     = sparseMatrix(c(1),j=c(ip4),x=c(-1),dims=c(1,np))
  Qp4sp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qp4)
  qp4sp   = cbind(Matrix(rep(0,n),nrow=1),qp4) 
  QC[[6]] = list(Qc=Qp4sp,q=toSparseVector(qp4sp),rhs=0,sense='=')


  # 
  Qp5     = sparseMatrix(ip3,ip2,x=1,dims=c(np,np))+
            sparseMatrix(ip4,ic,x=1,dims=c(np,np))+
            sparseMatrix(ip4,ip2,x=2,dims=c(np,np))
  qp5     = sparseMatrix(c(1),j=c(ip5),x=c(-1),dims=c(1,np))
  Qp5sp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qp5)
  qp5sp   = cbind(Matrix(rep(0,n),nrow=1),qp5) 
  QC[[7]] = list(Qc=Qp5sp,q=toSparseVector(qp5sp),rhs=0,sense='=')

  # 
  Qp6     = sparseMatrix(ip2,ip2,x=1,dims=c(np,np))+
            sparseMatrix(ib,ib,x=1,dims=c(np,np))
  qp6     = sparseMatrix(c(1),j=c(ip6),x=c(-1),dims=c(1,np))
  Qp6sp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qp6)
  qp6sp   = cbind(Matrix(rep(0,n),nrow=1),qp6) 
  QC[[8]] = list(Qc=Qp6sp,q=toSparseVector(qp6sp),rhs=0,sense='=')
  
  Qp7      = sparseMatrix(ip7,ip6,x=1,dims=c(np,np))
  Qp7sp    = bdiag(Matrix(rep(0,n*n),nrow=n),Qp7)
  QC[[9]] = list(Qc=Qp7sp,q=0,rhs=1,sense='=')

#  # optional bounds, hopefully help with performance
#  Qbnd     = sparseMatrix(ip7,ip2,x=1,dims=c(np,np))
#  Qbndsp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qbnd)
#  QC[[10]] = list(Qc=Qbndsp,q=0,rhs=1,sense='<')
#  QC[[11]] = list(Qc=Qbndsp,q=0,rhs=-1,sense='>')
#  
#  Qbnd     = sparseMatrix(ip7,ib,x=1,dims=c(np,np))
#  Qbndsp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qbnd)
#  QC[[12]] = list(Qc=Qbndsp,q=0,rhs=1,sense='<')
#  QC[[13]] = list(Qc=Qbndsp,q=0,rhs=-1,sense='>')
#  
#  Qbnd     = sparseMatrix(ib,ib,x=1,dims=c(np,np))
#  Qbndsp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qbnd)
#  QC[[14]] = list(Qc=Qbndsp,q=0,rhs=1e-5,sense='>')
  
  return(QC)
}