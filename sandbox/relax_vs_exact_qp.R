require(devtools)
require(gurobi)
load_all()
#runPcNCQ(freq=23.5,Nfine=144,Nmeas=48,con_mode='exact',WorkLimit=20,Threads=12)
Nmeas=48
Nfine=144
res=runPcNCQ(freq=23.5,Nfine=Nfine,Nmeas=Nmeas,con_mode='relax',WorkLimit=5,Threads=12,MIPFocus=2)
sol=res$sol
model=res$model


a=sol$x[Nfine+1]
b=sol$x[Nfine+2]
c=sol$x[Nfine+3]
p1=sol$x[Nfine+4]
p2=sol$x[Nfine+5]
p3=sol$x[Nfine+6]
p4=sol$x[Nfine+7]
p5=sol$x[Nfine+8]
p6=sol$x[Nfine+9]
p7=sol$x[Nfine+10]

print(p5*p7)
#expect_lte(p1,a^2-2*a*c+4*b^2+c^2)
#expect_gte(p2,a/2-c/2-p1/2)
#expect_gte(p3,a*p2)
#expect_lte(p5,p3*p2+p4*(c+2*p2))
#expect_gte(p6,p2^2+b^2)
#expect_lte(p7,1/p6)
print(matrix(c(a,b,b,c),nrow=2) |> eigen() |> (\(x) x$values)() |> min())

p1 = sqrt(a^2 - 2*a*c + 4*b^2 + c^2)
p2 = a/2 - c/2 - p1/2
p3 = a*p2 
p4 = b^2 
p5 = p3*p2 +p4*(c+2*p2) 
p6 = p2^2 + b^2
p7 = 1/p6

print(p5*p7)
print(sol$x[Nfine+8]<p5)
print(sol$x[Nfine+10]<p7)

#sol=runPcNCQ(freq=1,Nfine=80,Nmeas=36,con_mode='relax',check_start = T,WorkLimit=30,Threads=12)

