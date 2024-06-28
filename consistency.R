require(devtools)
require(dplyr)
load_all()

param=list(Amp=2,acro=runif(1)*2*pi,freq=1)
mt = c(0,1,2,4,5,6)/24
evalExactPower(mt,param)
evalMonteCarloPower(mt,param,1e5)

mt = c(1:24)/24 -1/24 
print(evalExactPower(mt,param))
print(evalMonteCarloPower(mt,param,1e5))
