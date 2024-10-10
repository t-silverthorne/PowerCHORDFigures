source('figures/fig_settings.R')
if (pub_qual){
  Nmc_glob=1e4
}else{
  Nmc_glob=1e3
}
motivation_fig=function(freq,equispaced_FN,Amp=2,Nm=24,Nm1=5,Nmc=Nmc_glob){
  tunif   = c(1:Nm)/Nm-1/Nm
  Nm2     = Nm-Nm1
  mt1     = c(1:Nm1)/Nm1-1/Nm1
  mt2     = c(1:Nm2)/Nm2-1/Nm2
  tirr    = c(mt1*.5, .5*mt2 +.5)
  acrovec = 2*pi*runif(Nmc)
  
  Yunif   = Amp*cos(outer(acrovec,2*pi*freq*tunif,'-'))+matrix(rnorm(length(tunif)*Nmc),nrow=Nmc)
  Yirr    = Amp*cos(outer(acrovec,2*pi*freq*tirr,'-'))+matrix(rnorm(length(tirr)*Nmc),nrow=Nmc)
  
  df_unif = data.frame(meas='equispaced',acro=acrovec,pval= matrixTests::row_cosinor(Yunif,tunif,1/freq) %>% {.$pvalue}) 
  df_irr  = data.frame(meas='irregular',acro=acrovec,pval= matrixTests::row_cosinor(Yirr,tirr,1/freq) %>% {.$pvalue}) 
  
  if (equispaced_FN){
    cmap_cust = c('true'=rgb(.05,0.5,.06),
                  'irregular'=rgb(.05,0.5,.06),
                  'equispaced'=rgb(.81,.06,.13))
    fn_ind  = sample(which(df_unif$pval>.05),1)
  }else{
    cmap_cust = c('true'=rgb(.05,0.5,.06),
                  'irregular'=rgb(.05,0.5,.06),
                  'equispaced'=rgb(.05,0.5,.06))
    fn_ind = sample(which(df_unif$pval<.05 & df_irr$pval<.05),1)
  }
  
  df_unif = df_unif[df_unif$pval<.05,]
  df_irr  = df_irr[df_irr$pval<.05,]
  
  df_all  = data.frame(meas='true',acro=acrovec,pval=NA)
  
  df=rbind(df_unif,df_irr,df_all)
  
  df$meas = factor(df$meas,levels=c('true','irregular','equispaced'))
  
  rad_brk = c(0,pi,2*pi)
  rad_lab = c(expression(0),
              expression(pi),
              expression(2*pi))
  plt= df %>% ggplot(aes(x=acro,fill=meas))+geom_histogram(aes(y=after_stat(density)))+facet_wrap(~meas)+
    scale_fill_manual(values=cmap_cust)
  
  plt = plt + scale_x_continuous(limits=c(0,2*pi),
                                 breaks =rad_brk,
                                 labels = rad_lab)
  plt = plt + theme(legend.position='none')
  plt = plt + labs(x=element_text('acrophase (rad)'),
                   y=element_text('density'))
  plt = plt + theme(
    strip.background=element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.25)
  )
  plt=plt+theme(text=element_text(size=fsize))
  pbot = plt
  
  
  tfine = seq(0,1,.001)
  
  acro = acrovec[fn_ind]
  dfine = data.frame(meas='true',time=tfine,signal=Amp*cos(2*pi*freq*tfine-acro))
  dunif = data.frame(meas='equispaced',time=tunif,signal=Yunif[fn_ind,])
  dirr  = data.frame(meas='irregular',time=tirr,signal=Yirr[fn_ind,])
  
  tdf = rbind(dfine,dunif,dirr)
  
  tdf$meas = factor(tdf$meas,levels=c('true','irregular','equispaced'))
  plt = tdf %>% ggplot(aes(x=time,y=signal,color=meas))+
    geom_line(data=tdf %>% filter(meas=='true'))+
    geom_point(data=tdf %>% filter(meas!='true'))+facet_wrap(~meas)+
    scale_color_manual(values=cmap_cust)+ theme(legend.position='none')
  plt = plt + labs(x=element_text('time'),
                   y=element_text('simulated signal'))
  plt = plt + scale_x_continuous(breaks=c(0,1))
  plt = plt + theme(
    strip.background=element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.25)
  )
  plt=plt+theme(text=element_text(size=fsize))
  ptop = plt
  
  Fig = ptop/pbot
  return(Fig) 
}

set.seed(28)

###########################
# Phase dependent power motivation 
###########################
p1=as.ggplot(motivation_fig(1,F,Nm = 8,Nm1=5))
p2=as.ggplot(motivation_fig(4,T,Nm=8,Nm1=5)) 


###########################
# Heatmap
###########################
if (pub_qual){
  Nfreq = 2^10
  Nacro = 2^7+1
}else{
  Nfreq = 2^9
  Nacro = 2^6+1
}
Amps = c(1)
Nmvec  = c(10,20,30)
acros=seq(0,2*pi,length.out=Nacro)
acros=acros[1:(length(acros)-1)]
freqs  = seq(1,24,length.out=Nfreq)
pars   = expand.grid(Amp=Amps,Nmeas=Nmvec,freq=freqs,acro=acros)

pars$power = c(1:dim(pars)[1]) %>% mclapply(mc.cores=8,function(ii){
  x=pars[ii,]
  Amp   = x$Amp
  Nmeas = x$Nmeas
  acro  = x$acro
  freq  = x$freq
  mt = c(1:Nmeas)/Nmeas-1/Nmeas
  return(evalExactPower(mt,param=list(Amp=Amp,acro=acro,freq=freq)))
})

pars$power = as.numeric(pars$power)

pars =pars %>% mutate(A=Amp,N=Nmeas)

rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),expression(pi/2),
            expression(pi),expression(3*pi/2),
            expression(2*pi))

p3 = pars %>% ggplot(aes(x=freq,y=acro,fill=power))+
  geom_raster()+
  facet_grid(A~N,labeller = purrr::partial(label_both, sep = " = "))+
  scale_y_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  scale_x_continuous(limits=c(0,24),
                                 breaks =seq(0,24,4),
                                 labels =seq(0,24,4))+
  scale_fill_viridis_c(limits=c(0,1))+
  labs(y='acrophase (rad)',x='frequency (cycles/day)')

plt_width=6
p3=p3+theme(legend.position='bottom',
              legend.key.width = unit(plt_width*.15, "in"),
              legend.title= element_text(hjust = 0.5),
              legend.direction = "horizontal")
p3 = p3 + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
p3=p3+theme(text=element_text(size=fsize))

Fig = p1/p2/p3 + plot_annotation(tag_levels='A')+plot_layout(heights=c(1.5,1.5,.75))
show_temp_plt(Fig,6,7)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'fig1.png'),
       Fig,
       width=6,height=7,
       device='png',
       dpi=600)
