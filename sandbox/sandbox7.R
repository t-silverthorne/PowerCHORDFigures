require(devtools)
load_all()
param = list(Amp=1, freq=1, acro=0)

mt = 1:24/24 - 1/24
print("24 pts covering full")
print(evalExactPower(mt,param))
print(evalMonteCarloPower(mt,param,1e6))


mt    = 0.5*(1:12/12 - 1/12)
print("12 pts covering first half")
print(evalExactPower(mt,param))
print(evalMonteCarloPower(mt,param,1e6))

