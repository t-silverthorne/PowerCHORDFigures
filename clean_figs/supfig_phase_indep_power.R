source('clean_figs/clean_theme.R')
Nfreq    = 2^10
Amp      = 1
mc_cores = 12


# only look at differential evolution results
#TODO: move to new dir
am = readRDS('clean_figs/data/powerCHORD_even_sols.RDS')
am = am[am@''$method=='diffEVCR',]
df = am@''

###########################
# power heatmap 
###########################
Nfreq_plt = 2^8
Nacro     = 2^6+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]
fmin      = .75 
fmax      = 1.25
freqs_plt = seq(fmin,fmax,length.out=Nfreq_plt)
pars=expand.grid(Amp=c(1),
                 Nmeas=c(24),
                 acro=acros,
                 freq=freqs_plt,
                 type=c('equispaced design'))

pwr_vec = c(1:dim(pars)[1]) %>% mclapply(mc.cores=12,function(ii){
  x      = pars[ii,]
  Nmeas  = as.numeric(x[['Nmeas']])
  acro   = as.numeric(x[['acro']])
  freq   = as.numeric(x[['freq']])
  Amp    = as.numeric(x[['Amp']]) 
  
  if(x$type == 'equispaced design'){
    mt = c(1:Nmeas)/Nmeas-1/Nmeas 
  }else if(x$type =='irregular design'){
    mt=as.numeric(am[am@''$Nmeas==Nmeas & am@''$fmin==fmin & am@''$fmax==fmax,])
    mt=mt[!is.na(mt)]
  }else{
    stop('unknown type')
  }
  if (length(mt)==Nmeas){
    power=evalPower(mt,Amp=Amp,freq=freq,acro=acro,design='equispaced')
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
  scale_y_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  labs(y='acrophase (rad)',x='frequency (cycles/day)')+
  guides(fill = guide_colorbar(barheight = unit(0.3, "inches")))+ 
  scale_fill_viridis_c(limits=c(0.6,.9),breaks=c(0.6,.7,.8,0.9))
plt=plt+clean_theme()
plt=plt+guides(fill=guide_colorbar(title.position='right'))
plt=plt+theme(legend.direction='vertical',
              legend.title = element_text(angle = 90,hjust=0.5))
phmap = plt
phmap

show_temp_plt(phmap,6,2)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'supfig_ipower.png'),
       phmap,
       width=6,height=2,
       device='png',
       dpi=600)
