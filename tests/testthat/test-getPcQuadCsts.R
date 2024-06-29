test_that("extra outputs match", {
  freq  = 0.59
  Nmeas = 23
  Nfine = 244
  tcvx  = runif(1)
  QC    = getPcQuadCsts(c(tcvx),freq,freq,1,Nfine,Nmeas) 

  expect_equal(QC[[1]]$q,c(rep(0,Nfine),-1))
  expect_equal(QC[[1]]$sense,'>')
})


test_that("agreement with one freq", {
  freq  =  100*runif(1) 
  Nmeas = sample(c(10:40),1)
  Nfine = sample((Nmeas+20):(Nmeas+50),1)
  tcvx  = runif(1)
  QC    = getPcQuadCsts(c(tcvx),freq,freq,1,Nfine,Nmeas) 

  # quadratic form evaluated on discrete grid using matrix from constraint
  tau = c(1:Nfine)/Nfine - 1/Nfine
  mu  = rep(0,Nfine)
  mu[sample(1:length(tau),Nmeas)]=1
  Binv = QC[[1]]$Qc
  Binv = Binv[c(1:Nfine),c(1:Nfine)]
  qform1 = t(mu)%*%Binv%*%mu

  # quadratic form evaluated directly from Fisher information matrix 
  mt = tau[mu>0]
  A  = matrix(c(0,0,1,0,0,1),nrow=3,byrow=T)
  X  = matrix(c(rep(1,Nmeas),cos(2*pi*freq*mt),sin(2*pi*freq*mt)),ncol=3)
  B  = t(A)%*%solve(t(X)%*%X,A)

  beta = matrix(c(sqrt(tcvx),sqrt(1-tcvx)),nrow=2)
  qform2 = t(beta)%*%solve(B,beta)
  expect_equal(qform1,qform2)
  
})

test_that("agreement with several freqs", {
  fmin  = 0.5
  fmax  = 24
  Nfreq = 30
  Nmeas = sample(c(10:40),1)
  Nfine = sample((Nmeas+20):(Nmeas+50),1)
  tcvxs = runif(Nfreq)
  QC    = getPcQuadCsts(tcvxs,fmin,fmax,Nfreq,Nfine,Nmeas) 

  # quadratic form evaluated on discrete grid using matrix from constraint
  tau = c(1:Nfine)/Nfine - 1/Nfine
  mu  = rep(0,Nfine)
  mu[sample(1:length(tau),Nmeas)]=1 # discrete time
  mt = tau[mu>0] # continuous time
  
  A  = matrix(c(0,0,1,0,0,1),nrow=3,byrow=T)
  qform1=c()
  qform2=c()
  fvec = seq(from=fmin,to=fmax,length.out=Nfreq)
  for(ii in c(1:Nfreq)){
    Binv = QC[[ii]]$Qc
    Binv = Binv[c(1:Nfine),c(1:Nfine)]
    qform1[ii] = t(mu)%*%Binv%*%mu

    tcvx = tcvxs[ii]
    freq = fvec[ii]
    X    = matrix(c(rep(1,Nmeas),cos(2*pi*freq*mt),sin(2*pi*freq*mt)),ncol=3)
    B    = t(A)%*%solve(t(X)%*%X,A)

    beta = matrix(c(sqrt(tcvx),sqrt(1-tcvx)),nrow=2)
    qform2[ii] = t(beta)%*%solve(B,beta)
  }
  expect_equal(qform1,qform2)
  
})