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
# Plot of power gain 
###########################
fdf=c(1:dim(am)[1]) %>% mclapply(mc.cores=12,function(ii){
  mt = as.numeric(am[ii,])
  mt = mt[!is.nan(mt)] 
  Nm = am@''[ii,]$Nmeas
  if (Nm == length(mt)){
     fmin=am@''[ii,]$fmin
     fmax=am@''[ii,]$fmax
     pwr = evalWorstPowerMultiFreq(mt,
                                   fmin=fmin,fmax=fmax,Nfreq=Nfreq,Amp=Amp,
                                   design='general')
     mt_unif = c(1:Nm)/Nm - 1/Nm 
     pwr_unif = evalWorstPowerMultiFreq(mt_unif,
                                        fmin=fmin,fmax=fmax,Nfreq=Nfreq,Amp=Amp,
                                        design='equispaced')
     return(cbind(am@''[ii,],
                  data.frame(pwr=pwr,pwr_unif=pwr_unif)))
  }else{
    stop('wrong number of measurement times')
  }
}) %>% rbindlist() %>% data.frame()

fdf$d_power = fdf$pwr-fdf$pwr_unif

# summary statistic
fdf$rel_gain = 100*(fdf$d_power/fdf$pwr_unif) |> abs()
fdf |> filter(d_power<0) |> select(rel_gain) |> max()

fdf = fdf |> mutate(dpower_plt = ifelse(abs(d_power)>.01,d_power,NA))
fdf = fdf %>% mutate(N=Nmeas)
plt = fdf %>% ggplot(aes(x=fmin,y=fmax,color=dpower_plt))+geom_point(size=3)+
  facet_wrap(~N,labeller = purrr::partial(label_both, sep = " = "),nrow=1)+
  scale_color_viridis_c(limits=c(0,1)) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12))+
  scale_y_continuous(breaks=c(4,8,12,16,20,24))+
  guides(fill = guide_colorbar(barheight = unit(0.4, "inches"))) 
plt=plt+clean_theme()
plt = plt + labs(x=element_text('minimum frequency (cycles/day)'),
                 y=element_text('maximum frequency\n(cycles/day)'),
                 color='power difference')
plt=plt+guides(color=guide_colorbar(title.position='right'))
plt=plt+theme(legend.direction='vertical',
                legend.title = element_text(angle = 90,hjust=0.5))
plt_gainWC = plt
plt


###########################
# power heatmap 
###########################
Nfreq_plt = 2^8
Nacro     = 2^6+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]
fmin      = 1
fmax      = 24
freqs_plt = seq(fmin,fmax,length.out=Nfreq_plt)
pars=expand.grid(Amp=c(1),
                 Nmeas=c(24,48),
                 acro=acros,
                 freq=freqs_plt,
                 type=c('equispaced design','irregular design'))

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
    if (x$type=='equispaced'){
      power=evalPower(mt,Amp=Amp,freq=freq,acro=acro)
    } else{
      power=evalPower(mt,Amp=Amp,freq=freq,acro=acro)
    }
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
  labs(y='acrophase (rad)',x='frequency (cycles/day)')+
  guides(fill = guide_colorbar(barheight = unit(0.3, "inches"))) 
plt=plt+clean_theme()
plt=plt+guides(fill=guide_colorbar(title.position='right'))
plt=plt+theme(legend.direction='vertical',
              legend.title = element_text(angle = 90,hjust=0.5))
phmap = plt
phmap

Fig1=plt_gainWC/phmap+ plot_annotation(tag_levels='A')+
  plot_layout(heights=c(1,1))

show_temp_plt(Fig1,6,4)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f2_broadprior1.png'),
       Fig1,
       width=6,height=3.5,
       device='png',
       dpi=600)
