require(devtools)
load_all()

freq  = 1
npar  = 8
Nfine = 144 
Nmeas = 24 

#require(gurobi)

tau   = c(1:Nfine)/Nfine-1/Nfine
model=list()
LC = getPcLinCsts(Nmeas,Nfine,Inf,npar=npar)

model$A          = LC$A
model$sense      = LC$sense_list %>% unlist() 
model$rhs        = LC$rhs_list %>% unlist() 
model$obj        = c(rep(0,Nfine+npar-1),1)
  
model$start      = getInitialState(freq=freq,Nfine=Nfine,Nmeas=Nmeas) 
model$modelsense = 'max'
model$modelname  = 'Power_CHORD'
model$vtype      = c(rep('B',Nfine),rep('C',npar))

Qmats            = getFourQuadBlocks(freq,Nfine,Nmeas)
model$quadcon    = getNonConvexQC(Nfine,
                   Matrix(Qmats$Cm11),
                   Matrix(Qmats$Cm12),
                   Matrix(Qmats$Cm22))

sol=gurobi(model,params=list(WorkLimit=20,Threads=12))