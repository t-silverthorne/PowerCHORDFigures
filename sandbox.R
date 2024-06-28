require(dplyr)
require(gurobi)

QQ = matrix(rnorm(n^2),n,n)
n = 50
Nmeas =4
model =list()
model$Q   = -QQ 
model$A   = matrix(rep(1,n),nrow=1)
model$rhs = Nmeas
model$vtype = rep('B',n)
model$modelsense = 'min'

gurobi(model)