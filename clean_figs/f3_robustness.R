source('clean_figs/clean_theme.R')
Nfreq    = 2^10
Amp      = 1
mc_cores = 12

# only look at differential evolution results
#am = readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
am = readRDS('clean_figs/data/powerCHORD_even_sols.RDS')
am = am[am@''$method=='diffEVCR',]
df = am@''

###########################
# Plot of raw solutions 
###########################
fmax=24
fmin=1
Nmeas_vals = c(24,48)

sloc     = am[am@''$Nmeas %in% Nmeas_vals & am@''$fmin==fmin & am@''$fmax==fmax,]
sloc[1,] = sloc[1,]-min(sloc[1,],na.rm=T)
sloc[2,] = sloc[2,]-min(sloc[2,],na.rm=T)
splt     = sloc %>% stack()
splt     = splt[!is.nan(splt$value),]
splt     = splt %>% mutate(N=Nmeas)

plt=splt[splt$Nmeas %in% Nmeas_vals, ] %>% ggplot(aes(x=value,y=0))+geom_point(size=.5)+
  facet_grid(N~.,
             labeller = purrr::partial(label_both, sep = " = "))
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(labels=seq(0,24,4),
  breaks=seq(0,1,4/24),
  limits=c(0,1))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
psol=plt
psol


###########################
# Robustness to measurement timing 
###########################
Nmvec  = c(24,32,48)
nrep   = 100 
scales = seq(1,30,2)/60/24
pars   = expand.grid(scale=scales,type=c('irregular','equispaced'),Nm=Nmvec)
sols   = readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
df=c(1:dim(pars)[1]) %>% mclapply(mc.cores=12,function(ii){
  sc   = pars[ii,]$scale
  type = pars[ii,]$type
  Nm   = pars[ii,]$Nm
  fmax_loc = Nm/2
  param  = list(Amp=1,fmin=1,fmax=fmax_loc,Nfreq=2^10)
  if (type=='equispaced'){
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

plt=df_grp %>% ggplot(aes(x=time,y=Power,group=type,color=type)) +geom_line() +
  geom_errorbar(aes(ymin=lower,ymax=upper))+
  facet_wrap(~Nm)
plt = plt+clean_theme()
plt = plt + labs(x=element_text('collection time standard deviation (minutes)'),
                 y=element_text('power'),
                 color='design')
plt = plt + theme(legend.position='bottom')
prob = plt

Fig=psol+prob + plot_layout(widths=c(2,3))+plot_annotation(tag_levels='A')

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f2_broadprior2.png'),
       Fig,
       width=6,height=2,
       device='png',
       dpi=600)
