require(gurobi)
# for each frequency: (non-convex MIQCP)
# mu'Af_{ij}*mu - sf_ij = 0
# [tf1 tf2]'[sf_11 sf_12 sf_21 sf_22][tf1 tf2] = hf
# max_{mu,t,s}(hf)
# Nfine + (2+2+1)*Nfreq variables
n = 1e2 
A1 = matrix(rnorm(n*n),nrow=n)
A2 = matrix(rnorm(n*n),nrow=n)
A3 = matrix(c(0,1,0,0,0,0,0,0,0),nrow=3)


A1 = rbind(cbind(A1,0,0,0),0,0,0)
A2 = rbind(cbind(A2,0,0,0),0,0,0)
A3 = rbind(matrix(rep(0,n*(n+3)),nrow=n),cbind(matrix(rep(0,n*3),nrow=3),A3))

model=list()
model$A     = matrix(c(rep(1,n),0,0,0),nrow=1)
model$sense = '='
model$rhs   = 5
model$quadcon = list(list(Qc=A1,q=c(rep(0,n),-1,0,0),rhs=0,sense='='),
                     list(Qc=A2,q=c(rep(0,n),0,-1,0),rhs=0,sense='='),
                     list(Qc=A3,q=c(rep(0,n),0,0,-1),rhs=0,sense='<')
)
model$obj = matrix(c(rep(0,n),0,0,1),nrow=1)
model$modelsense='min'
model$vtype  = c(rep('B',n),rep('C',3))

params=list(WorkLimit=100)
sol=gurobi(model,params)