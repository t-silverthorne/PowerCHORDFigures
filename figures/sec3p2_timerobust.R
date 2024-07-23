source('figures/fig_settings.R')
Nmvec  = c(24,32,48)
if (pub_qual){
  nrep   = 100 
  scales = seq(1,30,2)/60/24
}else{
  nrep   = 5
  scales = seq(1,30,2)/60/24
}
pars   = expand.grid(scale=scales,type=c('WCP','equispaced'),Nm=Nmvec)
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
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)

# might need
plt = plt + labs(x=element_text('collection time standard deviation'),
                 y=element_text('power'),
                 color='design')
plt = plt+theme(text=element_text(size=fsize))
Fig = plt

show_temp_plt(Fig,6,2)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'fig4.png'),
       Fig,
       width=6,height=2,
       device='png',
       dpi=600)