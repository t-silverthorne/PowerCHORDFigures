getNonConvexLC=function(Nmeas=48,Nfine=288,npar=1){
  if(Nmeas>Nfine){
    stop('Invalid constraint: Nmeas>Nfine, must increase Nfine')
  }
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

  ncstr = 3
  sense_list=list()
  rhs_list = list()
  A1=sparseMatrix(rep(1,Nfine),c(1:Nfine),
                x=rep(1,Nfine),dims=c(ncstr,Nfine+npar))
  A2=sparseMatrix(c(2,2,2,2),c(ia,ic,ip1,ip2),x=c(1/2,-1/2,-1/2,-1),dims=c(ncstr,npar))
  A2=cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),A2)
  A3 =cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),sparseMatrix(c(3),c(ip1),x=1,dims=c(ncstr,npar)))
  A =A1+A2+A3
  sense_list = list('=','=','>')
  rhs_list = list(Nmeas,0,0)

  return(list(A          = A,
              sense_list = sense_list,
              rhs_list   = rhs_list))
}