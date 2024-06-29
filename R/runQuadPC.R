runQuadPC=function(mu=48,vfreq=NULL,WLim=10,Nfine=144,fmin=1,fmax=24,Nfreq=49,drts=Inf,model=NULL,...){
  Nmeas = sum(mu)
  tau   = c(1:Nfine)/Nfine-1/Nfine
    
  # construct linear constraints 
  LC = getPcLinCsts(Nmeas,Nfine,drts)
  
  # construct quadratic constraints
  if (is.null(vfreq)){ # freq eigenvalues to define quad constraints
    vfreq=runif(Nfreq)
  }
  QC = getPcQuadCsts(vfreq,fmin,fmax,Nfreq,Nfine,Nmeas)
  
  # add linear constraints to gurobi model 
  model$A          = LC$A
  model$sense      = LC$sense_list %>% unlist() 
  model$rhs        = LC$rhs_list %>% unlist() 
  model$obj        = c(rep(0,Nfine),1)
  
  # starting vector 
  if (length(mu==1)){
    mu = rep(0,Nfine)
    mu[sample(1:Nfine,Nmeas,replace=F)]=1
  }

  # models ettings
  model$start      = c(mu,0) 
  model$modelsense = 'max'
  model$modelname  = 'Power_CHORD'
  model$vtype      = c(rep('B',Nfine),'C')
  
  model$quadcon=QC
  params = list(WorkLimit=10,...)
    
  return(gurobi(model,params))
}