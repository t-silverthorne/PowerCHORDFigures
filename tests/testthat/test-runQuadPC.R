require(devtools)
require(dplyr)
load_all()
Nfine   = 288
Niter   = 200
fmin    = 1 
fmax    = 24 
Nfreq   = 49 
Nftest  = 2^10 
WorkLim = 200 

vfreq = runif(Nfreq)
mu    = 40
tau   = c(1:Nfine)/Nfine-1/Nfine
pwr_now = rep(NaN,Niter)
travel  = rep(NaN,Niter)

for (iter in 1:Niter){
  qsys  = runQuadPC(mu=mu,vfreq=vfreq,WorkLim=WorkLim,Nfine=Nfine,
                    fmin=fmin,fmax=fmax,Nfreq=Nfreq,Threads=12)
  mu    = qsys$sol$x[1:Nfine]
  vfreq_old = vfreq
  vfreq = rep(NaN,Nfreq)
  for (ff in 1:Nfreq){ # update eigenfreqs
    Cm11  = qsys$Qmats[[ff]]$Cm11
    Cm12  = qsys$Qmats[[ff]]$Cm12
    Cm21  = qsys$Qmats[[ff]]$Cm21
    Cm22  = qsys$Qmats[[ff]]$Cm22

    a11   = t(mu)%*%Cm11%*%mu
    a12   = t(mu)%*%Cm12%*%mu
    a21   = t(mu)%*%Cm21%*%mu
    a22   = t(mu)%*%Cm22%*%mu

    eigs  = eigen(matrix(c(a11,a12,a21,a22),nrow=2))

    vfreq[ff]=abs(eigs$vectors[,which.min(eigs$values)][1])

  }
  mt = tau[mu>0]
  pwr_now[iter] = evalWorstPowerMultiFreq(mt,list(Amp=1,fmin=fmin,fmax=fmax,Nfreq=Nftest))
  travel[iter] =norm(matrix(vfreq-vfreq_old))
  cat(paste0('\n\n\nCURRENT POWER ',pwr_now[iter],'\n'))
  cat(paste0('Eigen travel: ',travel[iter],'\n\n\n'))
}
cat(paste0(pwr_now,'\n'))
cat(paste0(travel,'\n'))