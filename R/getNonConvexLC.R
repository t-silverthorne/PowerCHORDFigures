getNonConvexLC=function(Nmeas=48,Nfine=288,npar=10,con_mode='exact'){
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

  if (con_mode=='exact'){
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
  }else if(con_mode=='relax'){
    ncstr = 8
    sense_list=list()
    rhs_list = list()
    A1=sparseMatrix(rep(1,Nfine),c(1:Nfine),
                  x=rep(1,Nfine),dims=c(ncstr,Nfine+npar))
    A2=sparseMatrix(c(2,2,2,2),c(ia,ic,ip1,ip2),x=c(1/2,-1/2,-1/2,-1),dims=c(ncstr,npar))
    A2=cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),A2)
    A3 =cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),sparseMatrix(c(3),c(ip1),x=1,dims=c(ncstr,npar)))
    A4 =cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),sparseMatrix(c(4),c(ip2),x=1,dims=c(ncstr,npar)))
    A5 =cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),sparseMatrix(c(5),c(ip3),x=1,dims=c(ncstr,npar)))
    A6 =cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),sparseMatrix(c(6),c(ip5),x=1,dims=c(ncstr,npar)))
    A7 =cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),sparseMatrix(c(7),c(ip6),x=1,dims=c(ncstr,npar)))
    A8 =cbind(Matrix(rep(0,Nfine*ncstr),nrow=ncstr),sparseMatrix(c(8),c(ip7),x=1,dims=c(ncstr,npar)))
    A =A1+A2+A3+A4+A5+A6+A7+A8
    sense_list = list('=','<',
                     '>','<','<','>','>','>')
    rhs_list = list(Nmeas,0,0,0,0,0,0,0)

  }else{
    stop('unknown constraint mode')
  }
  return(list(A          = A,
              sense_list = sense_list,
              rhs_list   = rhs_list))
}