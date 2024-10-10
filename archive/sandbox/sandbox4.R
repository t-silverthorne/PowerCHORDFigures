require(devtools)
load_all()

freq  = 24
npar  = 10 
Nfine = 144 
Nmeas = 48 

#require(gurobi)

tau   = c(1:Nfine)/Nfine-1/Nfine
model=list()
LC = getNonConvexLC(Nmeas,Nfine,npar=npar)

model$A          = LC$A
model$sense      = LC$sense_list %>% unlist() 
model$rhs        = LC$rhs_list %>% unlist() 
model$obj        = rep(0,Nfine+npar)
  
model$start      = getInitialState(freq=freq,Nfine=Nfine,Nmeas=Nmeas) 
model$modelsense = 'max'
model$modelname  = 'Power_CHORD'
model$vtype      = c(rep('B',Nfine),rep('C',npar))

Qmats            = getFourQuadBlocks(freq,Nfine,Nmeas)
quadcon          = getNonConvexQC(Nfine,
                   Matrix(Qmats$Cm11),
                   Matrix(Qmats$Cm12),
                   Matrix(Qmats$Cm22))
model$quadcon    = quadcon
model$Q          = getNonConvexQobj(Nfine,npar)
model$start      = getInitialState(freq,Nfine,Nmeas)

check_start=F
if (check_start){
  model$lb         = model$start 
  model$ub         = model$start 
}else{
  model$lb         = rep(-1e1,Nfine+npar)
  model$ub         = rep(1e1,Nfine+npar) 
}

sol=gurobi(model,params=list(WorkLimit=1000,Threads=12))

tau = c(1:Nfine)/Nfine-1/Nfine
mt  = tau[sol$x[1:Nfine]>0]
cat(evalWorstPower(mt,list(freq=freq,Amp=1)),'\n')
cat(evalWorstPower(c(1:Nmeas)/Nmeas-1/Nmeas,list(freq=freq,Amp=1)),'\n')


#mu=sol$x[1:Nfine]
#
a   = sol$x[Nfine+1] # matrix element top
b   = sol$x[Nfine+2] # matrix element off-diag
c   = sol$x[Nfine+3] # matrix element bottom
p1  = sol$x[Nfine+4] # evec param 
p2  = sol$x[Nfine+5] # evec param 
p3  = sol$x[Nfine+6] # evec param 
p4  = sol$x[Nfine+7] # evec param 
p5  = sol$x[Nfine+8] # evec param 
p6  = sol$x[Nfine+9] # evec param 
p7  = sol$x[Nfine+10] # evec param 

ev  = eigen(matrix(c(a,b,b,c),nrow=2))


