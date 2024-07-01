require(devtools)
require(parallel)
require(data.table)
require(ggplot2)
load_all()

mode  = 2
Nfine = 288 
Nmeas = 30

if (mode==1){
  require(gurobi)
  sol   = runPcTrace(Nmeas,1,24,49*2,Nfine,WorkLimit=200,Threads=12)
  saveRDS(sol,'sandbox/temp.RDS')
}else{
  sol=readRDS('sandbox/temp.RDS')
  tau = c(1:Nfine)/Nfine - 1/Nfine

  mt      = tau[sol$x[1:Nfine]>0]
  mt_unif = c(1:Nmeas)/Nmeas-1/Nmeas

  Nfreq   = 2^6
  Nacro   = 2^6
  freqs   = seq(1,24,length.out=Nfreq)
  acro    = seq(1,2*pi,length.out=Nacro+1)
  acro    = acro[1:Nacro]

  pars=expand.grid(acro=acro,freq=freqs,type=c('unif','opt'))

  df=c(1:dim(pars)[1])%>% mclapply(mc.cores=12,function(ind){
    param=list(Amp=1,acro=pars[ind,]$acro,freq=pars[ind,]$freq)
    
    if (toString(pars[ind,]$type)=='unif'){
      pwr=evalExactPower(mt_unif,param)
    }else{
      pwr=evalExactPower(mt,param)
    }
    return(data.frame(pars[ind,],power=pwr))
  })%>%data.table::rbindlist() %>%data.frame()

  plt=df |> ggplot(aes(x=freq,y=acro,fill=power))+
    geom_raster()+
    scale_fill_viridis_c(limits=c(0,1))+
    facet_wrap(~type)
  plt
}
plt
print('done')