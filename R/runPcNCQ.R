runPcNCQ=function(freq=1,Nfine=288,Nmeas=48,con_mode='exact',check_start=F,...){
  npar=10
  tau   = c(1:Nfine)/Nfine-1/Nfine
  model = list()
  LC    = getNonConvexLC(Nmeas,Nfine,npar=npar,con_mode)

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
                    Matrix(Qmats$Cm22),con_mode=con_mode)
  model$quadcon    = quadcon
  model$Q          = getNonConvexQobj(Nfine,npar)
  model$start      = getInitialState(freq,Nfine,Nmeas)

  if (check_start){
    model$lb         = model$start 
    model$ub         = model$start 
  }else{
    model$lb         = rep(-Nmeas,Nfine+npar)
    model$ub         = rep(Nmeas,Nfine+npar) 
    #U=Nmeas
    #L=0.5
    #model$lb         = c(rep(0,Nfine),c(-U,-U,-U,L,-U,-U,L,L,L,L))
    #model$ub         = c(rep(1,Nfine),c( U, U, U,U,-L,-L,U,U,U,U))
  }

  sol=gurobi(model,params=list(...))
  return(list(sol=sol,model=model))
}
