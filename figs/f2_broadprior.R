source('figs/fig_settings.R')
Nfreq    = 2^10
Amp      = 1
mc_cores = 12
# aesthetics
run_aesthetics =T
plt_height = 3
plt_width  = 6

# only look at differential evolution results
am = readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
am = am[am@''$method=='diffEVCR',]
df = am@''


###########################
# Plot of power gain 
###########################
fdf=c(1:dim(am)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
  mt = as.numeric(am[ii,])
  mt = mt[!is.nan(mt)] 
  Nm = am@''[ii,]$Nmeas
  if (Nm == length(mt)){
     fmin=am@''[ii,]$fmin
     fmax=am@''[ii,]$fmax
     param=list(fmin=fmin,fmax=fmax,Nfreq=Nfreq,Amp=Amp)
     pwr = evalWorstPowerMultiFreq(mt,param)
     mt_unif = c(1:Nm)/Nm - 1/Nm 
     pwr_unif = evalWorstPowerMultiFreq(mt_unif,param)
     
     ncp = evalMinEigMultiFreq(mt,param)
     return(cbind(am@''[ii,],
                  data.frame(pwr=pwr,pwr_unif=pwr_unif,ncp=ncp)))
  }else{
    stop('wrong number of measurement times')
  }
}) %>% rbindlist() %>% data.frame()

fdf$d_power = fdf$pwr-fdf$pwr_unif

fdf = fdf %>% mutate(N=Nmeas)
plt = fdf %>% ggplot(aes(x=fmin,y=fmax,color=d_power))+geom_point(size=3)+
  facet_wrap(~N,labeller = purrr::partial(label_both, sep = " = "),nrow=1)+
  scale_color_viridis_c(limits=c(.01,1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12))+
  scale_y_continuous(breaks=c(4,8,12,16,20,24))
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt = plt + labs(x=element_text('minimum frequency'),
                 y=element_text('maximum frequency'),
                 color='power difference')
plt=plt+guides(color=guide_colorbar(title.position='right'))
plt=plt+theme(text=element_text(size=fsize),legend.direction='vertical',
                legend.title = element_text(angle = 90))
plt_gainWC = plt
plt

###########################
# Plot of raw solutions 
###########################
fmax=24
fmin=1
Nmeas_vals = c(24,48)

sloc = am[am@''$Nmeas %in% Nmeas_vals & am@''$fmin==fmin & am@''$fmax==fmax,]
sloc %>% stack()
splt = sloc %>% stack()
splt = splt[!is.nan(splt$value),]
splt
splt = splt %>% mutate(N=Nmeas)

plt=splt[splt$Nmeas %in% Nmeas_vals, ] %>% ggplot(aes(x=value,y=0))+geom_point(size=.5)+
  facet_grid(N~.,
             labeller = purrr::partial(label_both, sep = " = "))
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (days)'))
plt = plt+scale_x_continuous(limits=c(0,1),breaks=c(0,1))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+theme(text=element_text(size=fsize))
psol=plt

###########################
# power heatmap 
###########################
rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),
            expression(pi/2),
            expression(pi),
            expression(3*pi/2),
            expression(2*pi))
Nfreq_plt = 2^8
Nacro     = 2^6+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]
freqs_plt = seq(1,24,length.out=Nfreq_plt)
pars=expand.grid(Amp=c(1),
                 Nmeas=Nmeas_vals,
                 acro=acros,
                 freq=freqs_plt,
                 type=c('equispaced design','WCP design'))

pwr_vec = c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
  # unpack
  x      = pars[ii,]
  Nmeas  = as.numeric(x[['Nmeas']])
  acro   = as.numeric(x[['acro']])
  freq   = as.numeric(x[['freq']])
  Amp    = as.numeric(x[['Amp']]) 
  param =list(Amp=Amp,freq=freq,acro=acro) 
  
  if(x$type == 'equispaced design'){
    mt = c(1:Nmeas)/Nmeas-1/Nmeas 
  }else if(x$type =='WCP design'){
    mt=as.numeric(am[am@''$Nmeas==Nmeas & am@''$fmin==fmin & am@''$fmax==fmax,])
    mt=mt[!is.na(mt)]
  }else{
    stop('unknown type')
  }
  if (length(mt)==Nmeas){
    power=evalExactPower(mt,param)
  }else{
    stop('wrong length for measurement vector')
  }
  return(power)
})
pars$power = pwr_vec

pars$freq =as.numeric(pars$freq)
pars$acro =as.numeric(pars$acro)
pars$power=as.numeric(pars$power)

pars$N = factor(pars$Nmeas,unique(pars$Nmeas),paste0('N = ',unique(pars$Nmeas)))
plt = pars%>%
  ggplot(aes(x=freq,y=acro,fill=power))+
  geom_raster()+
  facet_grid(N~type)+
  scale_y_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  scale_fill_viridis_c(limits=c(0,1))+
  labs(y='acrophase (rad)',x='frequency (cycles/day)')

#plt=plt+guides(fill=guide_colorbar(title.position='left'))
#plt=plt+theme(legend.position='bottom',
#              legend.key.width = unit(plt_width*.15, "in"),
#              legend.direction = "horizontal")
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+guides(fill=guide_colorbar(title.position='right'))
plt=plt+theme(text=element_text(size=fsize),legend.direction='vertical',
                legend.title = element_text(angle = 90,hjust=0.5))
plt=plt+theme(text=element_text(size=fsize))
phmap = plt
phmap


###########################
# Robustness to measurement timing 
###########################
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
plt = plt + theme(legend.position='bottom')
prob = plt

Fig=plt_gainWC/psol/phmap/prob + plot_annotation(tag_levels='A')+
  plot_layout(heights=c(1.5,.8,1,1))

show_temp_plt(Fig,6,7)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f2_broadprior.png'),
       Fig,
       width=6,height=7,
       device='png',
       dpi=600)