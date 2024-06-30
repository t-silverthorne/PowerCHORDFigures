getNonConvexQobj=function(Nfine,npar){
  ip5 = 8
  ip7 = 10
  return(bdiag(Matrix(rep(0,Nfine*Nfine),nrow=Nfine),
         sparseMatrix(ip5,ip7,x=1,dims=c(npar,npar))))
}