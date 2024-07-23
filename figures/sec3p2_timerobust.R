require(dplyr)
require(data.table)
require(ggplot2)
require(ggplotify)
require(patchwork)
require(devtools)
require(parallel)
require(annmatrix)
load_all()
Nmvec  = c(24,32,48)
nrep   = 10  
scales = c(1:30)/60/24

pars   = expand.grid(scale=scales,type=c('WCP','uniform'),Nm=Nmvec)
param  = list(Amp=1,fmin=1,fmax=24,Nfreq=2^10)
sols   = readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
Nm=32
df=c(1:dim(pars)[1]) %>% mclapply(mc.cores=8,function(ii){
  sc   = pars[ii,]$scale
  type = pars[ii,]$type
  Nm   = pars[ii,]$Nm
  fmax_loc = Nm/2
  if (type=='uniform'){
    mt = c(1:Nm)/Nm -1/Nm 
  }else{
    filt = sols@Nmeas==Nm & sols@fmin == param$fmin & sols@fmax==fmax_loc & sols@method=='diffEVCR'
    mt = sols[filt,]
    mt = as.numeric(mt)
    mt = mt[!is.nan(mt)]
    if(length(mt)!=Nm){
      stop('wrong length meas vec')
    }
  }
  cbind(pars[ii,],
        data.frame(power=replicate(nrep,{evalWorstPowerMultiFreq(mt+rnorm(Nm,0,sd=sc),param=param)})))
  }
) %>% rbindlist() %>% data.frame()

df_grp = df %>% group_by(scale,type,Nm) %>% summarize(Power=mean(power),
                                     lower=quantile(power)[2],
                                     upper=quantile(power)[4])
df_grp$time = df_grp$scale*24*60

df_grp %>% ggplot(aes(x=time,y=Power,group=type,color=type)) +geom_line() +
  geom_errorbar(aes(ymin=lower,ymax=upper))+
  facet_wrap(~Nm)
