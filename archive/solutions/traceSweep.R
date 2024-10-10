require(devtools)
require(gurobi)
load_all()
fmins     = c(1,seq(2,12,2))
fmaxs     = c(1,seq(2,24,2))
Nmeas     = c(16,24,32,48)

pars      = expand.grid(fmin=fmins,fmax=fmaxs,Nmeas=Nmeas)
pars      = pars[pars$fmin<=pars$fmax,]

Nfine     = 144
Nfreq     = 49
WorkLimit = 60
Threads   = 12

sols_master = c(1:dim(pars)[1])%>%lapply(function(ind){
  fmin  = pars[ind,]$fmin
  fmax  = pars[ind,]$fmax
  Nmeas = pars[ind,]$Nmeas
  sol   = runPcTrace(Nmeas=Nmeas,fmin=fmin,fmax=fmax,Nfreq=Nfreq,Nfine=Nfine,
              WorkLimit=WorkLimit,Threads=Threads)
  
  return(list(pars=pars[ind,],sol=sol))
})

saveRDS(sols_master,'solutions/traceSweep_quick.RDS')