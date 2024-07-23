source('figures/fig_settings.R')
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
  scale_color_viridis_c(limits=c(.01,1)) 
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
                legend.title = element_text(angle = 90)
              )
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

plt=plt+guides(fill=guide_colorbar(title.position='left'))
plt=plt+theme(legend.position='bottom',
              legend.key.width = unit(plt_width*.15, "in"),
              legend.direction = "horizontal")
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+theme(text=element_text(size=fsize))
phmap = plt
phmap

Fig=plt_gainWC/psol/phmap + plot_annotation(tag_levels='A')

show_temp_plt(Fig,6,6)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'fig3.png'),
       Fig,
       width=6,height=6,
       device='png',
       dpi=600)
# aesthetics
