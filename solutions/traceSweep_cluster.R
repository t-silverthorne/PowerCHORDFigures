require(devtools)
require(parallel)
require(gurobi)
load_all()
fmins     = c(1,seq(2,12,1))
fmaxs     = c(1,seq(2,24,1))
Nmeas     = c(16,24,32,48)

pars      = expand.grid(fmin=fmins,fmax=fmaxs,Nmeas=Nmeas)
pars      = pars[pars$fmin<=pars$fmax,]

Nfine     = 144
Nfreq     = 49
WorkLimit = 60
Threads   = 1

sols_master = c(1:dim(pars)[1])%>%mclapply(
  mc.cores=as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK")),
function(ind){
  fmin  = pars[ind,]$fmin
  fmax  = pars[ind,]$fmax
  Nmeas = pars[ind,]$Nmeas
  sol   = runPcTrace(Nmeas=Nmeas,fmin=fmin,fmax=fmax,Nfreq=Nfreq,Nfine=Nfine,
              WorkLimit=WorkLimit,Threads=1)
  return(list(pars=pars[ind,],sol=sol))
})

saveRDS(sols_master,'solutions/traceSweep_hires.RDS')