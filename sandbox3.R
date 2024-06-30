require(Matrix)
n   = 288

# dummy q-form matrices, will be replaced
A11 = Matrix::.dense2sparse(Matrix(rnorm(n*n),nrow=n))
A12 = Matrix::.dense2sparse(Matrix(rnorm(n*n),nrow=n))
A22 = Matrix::.dense2sparse(Matrix(rnorm(n*n),nrow=n))

# indices for constraint forms
ia   = 1 # matrix element top
ib   = 2 # matrix element off-diag
ic   = 3 # matrix element bottom
ip1  = 4 # first evec nonlin
ip2  = 5 # second evec nonlin
iv   = 6 # eigenvec
ivs  = 7 # for square of v 
imx  = 8 # for full quadratic form
npar = 8
QC=list()

# construct the discrete quadratic forms
A11sp = bdiag(A11,Matrix(rep(0,npar*npar),nrow=npar))
A12sp = bdiag(A12,Matrix(rep(0,npar*npar),nrow=npar))
A22sp = bdiag(A22,Matrix(rep(0,npar*npar),nrow=npar))

q11 = sparseMatrix(1,j=ia,x=-1,dims=c(1,npar))
q12 = sparseMatrix(1,j=ib,x=-1,dims=c(1,npar))
q22 = sparseMatrix(1,j=ic,x=-1,dims=c(1,npar))

q11sp = cbind(Matrix(rep(0,n),nrow=1),q11) 
q12sp = cbind(Matrix(rep(0,n),nrow=1),q12) 
q22sp = cbind(Matrix(rep(0,n),nrow=1),q22) 

QC[[1]]=list(Qc=A11sp,q=q11sp,rhs=0,sense='=')
QC[[2]]=list(Qc=A12sp,q=q12sp,rhs=0,sense='=')
QC[[3]]=list(Qc=A22sp,q=q22sp,rhs=0,sense='=')

# parameterize first part of eigenvector
Qp1     = sparseMatrix(ia,ia,x=1,dims=c(npar,npar))+
          sparseMatrix(ia,ic,x=-2,dims=c(npar,npar))+
          sparseMatrix(ib,ib,x=4,dims=c(npar,npar))+
          sparseMatrix(ic,ic,x=-1,dims=c(npar,npar))+
          sparseMatrix(ip1,ip1,x=-1,dims=c(npar,npar))
Qp1sp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qp1)
QC[[4]] = list(Qc=Qp1sp,q=0,rhs=0,sense='=')

# parameterize second part of eigenvector
Qp2     = sparseMatrix(ib,ip2,x=1,dims=c(npar,npar))
Qp2sp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qp2)
QC[[5]] = list(Qc=Qp2sp,q=0,rhs=1,sense='=')

#TODO: check sign
# parameterize eigenvector
Qe      = sparseMatrix(ip1,ip2,x=-1/2,dims=c(npar,npar))+
          sparseMatrix(ic,ip2,x=-1,dims=c(npar,npar))
qe      = sparseMatrix(1,iv,x=-1,dims=c(1,npar))+
          sparseMatrix(1,ia,x=1/2,dims=c(1,npar))+
          sparseMatrix(1,ic,x=1/2,dims=c(1,npar))
Qesp    = bdiag(Matrix(rep(0,n*n),nrow=n),Qe)
qesp    = cbind(Matrix(rep(0,n),nrow=1),qe) 
QC[[6]] = list(QC=Qesp,q=qesp,rhs=0,sense='=')


# slack variable for square of eigenvector component 1
Qsl     = sparseMatrix(iv,iv,x=1,dims=c(npar,npar))
qsl     = sparseMatrix(1,ivs,x=-1,dims=c(1,npar))
Qslsp   = bdiag(Matrix(rep(0,n*n),nrow=n),Qsl)
qslsp   = cbind(Matrix(rep(0,n),nrow=1),qsl) 
QC[[7]] = list(QC=Qslsp,q=qslsp,rhs=0,sense='=')

# build full quadratic form
Qful    = sparseMatrix(ia,ivs,x=1,dims=c(npar,npar))+
          sparseMatrix(iv,ib,x=2,dims=c(npar,npar))
qful    = sparseMatrix(1,ic,x=1,dims=c(1,npar))+
          sparseMatrix(1,imx,x=-1,dims=c(1,npar))
Qfulsp  = bdiag(Matrix(rep(0,n*n),nrow=n),Qful)
qfulsp  = cbind(Matrix(rep(0,n),nrow=1),qful) 
QC[[8]] = list(QC=Qfulsp,q=qfulsp,rhs=0,sense='>') #TODO check direction

print(QC)

