source('PLOSfigures/clean_theme.R')
Nfreq    = 2^10
Amp      = 1


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

pwr_vec = c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
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
    power=evalPower(mt,Amp=Amp,freq=freq,acro=acro)
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


###########################
# near-Nyquist power fluctuations 
###########################

Nmin  = 8  
Nmax  = 40 
Nvals = Nmin:Nmax
Nfr   = 5e2

df=Nvals |> lapply(function(N){
  mt   = c(1:N)/N-1/N
  data.frame(N=N,
             freq=seq(N*0.9/2,N*1.1/2,length.out=Nfr)*2/N,
             power=evalWorstPowerMultiFreq(mt,fmin=0.9*N/2,fmax=1.1*N/2,
                                           Nfreq = Nfr,
                                           Amp=1,
                                           returnType='all')
  )
}) |> rbindlist() |> data.frame()


lplt=df |> ggplot(aes(x=freq,y=power,color=N,group=N))+geom_line()+
  scale_color_viridis_c(option='A')+clean_theme()+
  labs(y='power',x='relative frequency (f/Nyquist)',
       color='sample size')
lplt=lplt+guides(color=guide_colorbar(title.position='right'))
lplt=lplt+theme(legend.direction='vertical',
              legend.title = element_text(angle = 90,hjust=0.5))


Fig = (phmap|lplt)+plot_annotation(tag_levels='A')

ggsave('PLOSfigures/suppfig2.png',
       Fig,
       width=6,height=2,
       device='png',
       dpi=600)
