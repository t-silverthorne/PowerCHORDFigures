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
# Plot of power gain 
###########################
fdf=c(1:dim(am)[1]) %>% mclapply(mc.cores=8,function(ii){
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
     return(cbind(am@''[ii,],
                  data.frame(pwr=pwr,pwr_unif=pwr_unif)))
  }else{
    stop('wrong number of measurement times')
  }
}) %>% rbindlist() %>% data.frame()

fdf$d_power = fdf$pwr-fdf$pwr_unif

fdf = fdf |> mutate(dpower_plt = ifelse(abs(d_power)>.01,d_power,NA))
fdf = fdf %>% mutate(N=Nmeas)
plt = fdf %>% ggplot(aes(x=fmin,y=fmax,color=dpower_plt))+geom_point(size=3)+
  facet_wrap(~N,labeller = purrr::partial(label_both, sep = " = "),nrow=1)+
  scale_color_viridis_c(limits=c(.01,1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12))+
  scale_y_continuous(breaks=c(4,8,12,16,20,24))
plt=plt+clean_theme()
plt = plt + labs(x=element_text('minimum frequency'),
                 y=element_text('maximum frequency'),
                 color='power difference')
plt=plt+guides(color=guide_colorbar(title.position='right'))
plt=plt+theme(legend.direction='vertical',
                legend.title = element_text(angle = 90))
plt_gainWC = plt
plt


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
                 Nmeas=c(24,48),
                 acro=acros,
                 freq=freqs_plt,
                 type=c('equispaced design','irregular design'))

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
  }else if(x$type =='irregular design'){
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
plt = pars |> 
  ggplot(aes(x=freq,y=acro,fill=power))+
  geom_raster()+
  facet_grid(N~type)+
  scale_y_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  scale_fill_viridis_c(limits=c(0,1))+
  labs(y='acrophase (rad)',x='frequency (cycles/day)')
plt=plt+clean_theme()
plt=plt+guides(fill=guide_colorbar(title.position='right'))
plt=plt+theme(legend.direction='vertical',
              legend.title = element_text(angle = 90,hjust=0.5))
phmap = plt
phmap

plt_gainWC = plt_gainWC + theme(
  legend.key.height= unit(dev.size()[2]/25, "inches"))

phmap =  phmap + theme(
  legend.key.height= unit(dev.size()[2]/25, "inches"))

Fig1=plt_gainWC/phmap+ plot_annotation(tag_levels='A')+
  plot_layout(heights=c(1.5,1))

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f2_broadprior1.png'),
       Fig1,
       width=6,height=3,
       device='png',
       dpi=600)
