require(devtools)
require(gurobi)
require(dplyr)
require(parallel)

load_all()
Nmeas = 48 
Nfreq = 49
Nfine = 144
drts  = Inf
fmin  = 1
fmax  = 24
Niter = 3

c(1:10) %>% mclapply(mc.cores=10,function(ind){
  
  for (jj in c(1:Niter)){
  xfreq = runif(Nfreq)
  
  model=list()
  
  # construct linear constraints 
  LC = getPcLinCsts(Nmeas,Nfine,drts)
  
  # construct quadratic constraints
  QC = getPcQuadCsts(xfreq,fmin,fmax,Nfreq,Nfine)
  
  # add linear constraints to gurobi model 
  model$A          = LC$A
  model$sense      = LC$sense_list %>% unlist() 
  model$rhs        = LC$rhs_list %>% unlist() 
  model$obj        = c(rep(0,Nfine),1)
  
  # model settings
  mu_start = rep(0,Nfine)
  mu_start[sample(c(1:Nfine),Nmeas)]=1
  model$start      = c(mu_start,0) 
  model$modelsense = 'max'
  model$modelname  = 'Power_CHORD'
  model$vtype      = c(rep('B',Nfine),'C')
  
  model$quadcon=QC
  
  # gurobi settings
  params = list(Presolve=2,MIPFocus=3,NumericFocus=3,Threads=1,WorkLimit=10)
  
  # solve
  sol=gurobi(model,params)
    
  } 
})
