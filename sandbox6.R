
freq  = 4*runif(1)
Nfine = sample(100:200,1)
Nmeas = sample(10:40,1) 
s0=getInitialState(freq,Nfine,Nmeas)

a  = s0[Nfine+1]
b  = s0[Nfine+2]
c  = s0[Nfine+3]
p1 = s0[Nfine+4]
p2 = s0[Nfine+5]
p3 = s0[Nfine+6]
p4 = s0[Nfine+7]
p5 = s0[Nfine+8]
p6 = s0[Nfine+9]
p7 = s0[Nfine+10]
print(c+2*p2)
