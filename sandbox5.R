require(gurobi)
load_all()
fmin        = 1 
fmax        = 12 
Nfreq       = 24 
Nfine       = 288 
Nmeas       = 32 
npar        = 1
fvec        = seq(fmin,fmax,length.out=Nfreq)
model       = list()
model$A     = matrix(c(rep(1,Nfine),rep(0,npar)),nrow=1)
model$sense = '='
model$rhs   = Nmeas

model$obj        = c(rep(0,Nfine),1)
model$quadcon    = list()

for (ii in c(1:Nfreq)){
  Cm=getFourQuadBlocks(fvec[ii],Nfine,Nmeas)
  model$quadcon[[ii]]    = list(
                          Qc=bdiag(Cm$Cm11+Cm$Cm22,matrix(rep(-1,npar*npar),nrow=npar)),
                          rhs=0,
                          sense='>')

}
model$modelsense = 'max'
model$modelname  = 'Power_CHORD'
model$vtype      = c(rep('B',Nfine),rep('C',npar))

sol = gurobi(model,list(WorkLimit=1,Threads=12))
tau = c(1:Nfine)/Nfine - 1/Nfine
cat('Optimal: ',evalWorstPowerMultiFreq(tau[sol$x[1:Nfine]>0],
  param=list(Amp=1,fmin=fmin,fmax=fmax,Nfreq=2^10)),'\n')
cat('Uniform: ',evalWorstPowerMultiFreq(c(1:Nmeas)/Nmeas-1/Nmeas,
  param=list(Amp=1,fmin=fmin,fmax=fmax,Nfreq=2^10)),'\n')


print(as.numeric(sol$x[1:Nfine]>0))
