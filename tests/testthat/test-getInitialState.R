require(Matrix)
test_that("initial state satisfies Qforms", {
  freq  = runif(1)*3 
  Nfine = sample(100:200,1) 
  Nmeas = sample(20:40,1)
  npar  = 10 

  s0   = getInitialState(freq,Nfine,Nmeas)
  Qmat = getFourQuadBlocks(freq,Nfine,Nmeas)
  Q11  = Matrix(Qmat$Cm11)
  Q12  = Matrix(Qmat$Cm12)
  Q22  = Matrix(Qmat$Cm22)
  QC = getNonConvexQC(Nfine,Q11,Q12,Q22)
  
  for (ii in c(1:length(QC))){
    qloc    = QC[[ii]]
    if(length(qloc$Qc)==1){
      Qmat = Matrix(matrix(rep(0,(Nfine+npar)^2),nrow=Nfine+npar)) 
    }else{
      Qmat = qloc$Qc
    }
    if(length(qloc$q)==1){
      qlin = Matrix(rep(qloc$q,Nfine+npar),nrow=Nfine+npar)
    }else{
      qlin=Matrix(qloc$q)
    }
    lhs_num = as.numeric(Matrix::t(s0)%*%Qmat%*%s0+Matrix::t(qlin)%*%s0)
    rhs_num = as.numeric(qloc$rhs)
    
    if (qloc$sense == '='){
      expect_equal(lhs_num,rhs_num)
    }else if(qloc$sense=='>'){
      expect_gte(lhs_num,rhs_num)
    }
  }

})

test_that('state correctly encodes eigenvector',{
  freq  = runif(1)*10
  Nfine = sample(100:200,1) 
  Nmeas = sample(20:40,1)
  npar  = 10 

  s0   = getInitialState(freq,Nfine,Nmeas)
  Qmat = getFourQuadBlocks(freq,Nfine,Nmeas)
  Q11  = Matrix(Qmat$Cm11)
  Q12  = Matrix(Qmat$Cm12)
  Q22  = Matrix(Qmat$Cm22)
  QC = getNonConvexQC(Nfine,Q11,Q12,Q22)

  mu = s0[1:Nfine]

  a11=as.numeric(t(mu)%*%Q11%*%mu)
  a12=as.numeric(t(mu)%*%Q12%*%mu)
  a22=as.numeric(t(mu)%*%Q22%*%mu)
  A  = matrix(c(a11,a12,a12,a22),nrow=2)
  lhs=eigen(matrix(c(a11,a12,a12,a22),nrow=2)) %>% {.$values} %>% min()
  rhs=s0[Nfine+npar-2]*s0[Nfine+npar]

  Qobj =getNonConvexQobj(Nfine,npar)
  expect_equal(lhs,rhs)
})

test_that('sign condition on components of initial state',{
  freq  = runif(1)*20
  Nfine = sample(100:200,1) 
  Nmeas = sample(20:40,1)
  npar  = 10 

  s0   = getInitialState(freq,Nfine,Nmeas)

  s0[1:Nfine]
  expect_gte(s0[Nfine+1],0)
  #expect_lte(s0[Nfine+2],0) # INDEFINITE
  expect_gte(s0[Nfine+3],0)
  expect_gte(s0[Nfine+4],0)
  expect_lte(s0[Nfine+5],0) # can prove, factor square
  expect_gte(s0[Nfine+7],0)
  expect_gte(s0[Nfine+8],0)
  expect_gte(s0[Nfine+9],0)
  expect_gte(s0[Nfine+10],0)
})