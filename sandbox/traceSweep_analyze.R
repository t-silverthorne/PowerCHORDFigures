require(devtools)
require(data.table)
require(ggplot2)
require(parallel)
require(patchwork)
sols  = readRDS('solutions/traceSweep_hires.RDS')
load_all()
Nfine = 144
tau   = c(1:Nfine)/Nfine-1/Nfine


df = c(1:length(sols)) |> mclapply(mc.cores=12,function(ind){
  sol      = sols[[ind]]
  pars     = sols[[ind]]$pars

  mt_opt   = tau[sol$sol$x[1:Nfine]>0]
  mt_unif  = c(1:pars$Nmeas)/pars$Nmeas-1/pars$Nmeas 

  pwr_opt  = evalWorstPowerMultiFreq(mt_opt,param=list(Amp=1,
                fmin=pars$fmin,fmax=pars$fmax*.9999,Nfreq=2^6))
  pwr_unif = evalWorstPowerMultiFreq(mt_unif,param=list(Amp=1,
                fmin=pars$fmin,fmax=pars$fmax*.9999,Nfreq=2^6))
  gain     = pwr_opt-pwr_unif
  
  cbind(pars,data.frame(
    pwr_unif = pwr_unif,
    pwr_opt  = pwr_opt,
    gain     = gain,
    runtime  = sol$sol$runtime,
    status   = sol$sol$status,
    mipgap   = sol$sol$mipgap
  ))
}) |> rbindlist() |> data.frame()


p1=df |>  ggplot(aes(x=fmin,y=fmax,color=gain))+geom_point(size=3)+
    facet_wrap(~Nmeas,nrow=1)+
  scale_color_viridis_c(limits=c(0,1)) 

p2=df |>  ggplot(aes(x=fmin,y=fmax,color=pwr_opt))+geom_point(size=3)+
    facet_wrap(~Nmeas,nrow=1)+
  scale_color_viridis_c(limits=c(0,1)) 

p1/p2 + plot_annotation(tag_levels='A')

df[df$Nmeas==48 & df$fmax==24,]
