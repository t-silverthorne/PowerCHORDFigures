source('figs/fig_settings.R')
pub_qual=T
if (pub_qual){
  Nmc_glob=1e4
}else{
  Nmc_glob=1e3
}
require(matrixTests)
require(ggplot2)
require(data.table)
require(tidyr)
require(patchwork)
# load in cutsdp solutions
n=48
tau = c(1:n)/n -1/n
Xraw=read.csv2('matlab/LMI_formulation/cutsdp_sols.csv',header = F,sep=',')
head(Xraw)

motivation_fig=function(freq,equispaced_FN,Amp=2,Nm=24,Nm1=5,Nmc=Nmc_glob){
  tunif   = c(1:Nm)/Nm-1/Nm
  tirr    = tau[Xraw[3,]>1e-12]
  acrovec = 2*pi*runif(Nmc)
  
  Yunif   = Amp*cos(outer(acrovec,2*pi*freq*tunif,'-'))+matrix(rnorm(length(tunif)*Nmc),nrow=Nmc)
  Yirr    = Amp*cos(outer(acrovec,2*pi*freq*tirr,'-'))+matrix(rnorm(length(tirr)*Nmc),nrow=Nmc)
  
  df_unif = data.frame(meas='equispaced',acro=acrovec,pval= matrixTests::row_cosinor(Yunif,tunif,1/freq) %>% {.$pvalue}) 
  df_irr  = data.frame(meas='optimal',acro=acrovec,pval= matrixTests::row_cosinor(Yirr,tirr,1/freq) %>% {.$pvalue}) 
  
  if (equispaced_FN){
    cmap_cust = c('true'=rgb(.05,0.5,.06),
                  'optimal'=rgb(.05,0.5,.06),
                  'equispaced'=rgb(.81,.06,.13))
    fn_ind  = sample(which(df_unif$pval>.05),1)
  }else{
    cmap_cust = c('true'=rgb(.05,0.5,.06),
                  'optimal'=rgb(.05,0.5,.06),
                  'equispaced'=rgb(.05,0.5,.06))
    fn_ind = sample(which(df_unif$pval<.05 & df_irr$pval<.05),1)
  }
  
  df_unif = df_unif[df_unif$pval<.05,]
  df_irr  = df_irr[df_irr$pval<.05,]
  
  df_all  = data.frame(meas='true',acro=acrovec,pval=NA)
  
  df=rbind(df_unif,df_irr,df_all)
  
  df$meas = factor(df$meas,levels=c('true','equispaced','optimal'))
  
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
  plt = plt + labs(x=element_text('acrophase'),
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
  dirr  = data.frame(meas='optimal',time=tirr,signal=Yirr[fn_ind,])
  
  tdf = rbind(dfine,dunif,dirr)
  
  tdf$meas = factor(tdf$meas,levels=c('true','equispaced','optimal'))
  plt = tdf %>% ggplot(aes(x=time,y=signal,color=meas))+
    geom_line(data=tdf %>% filter(meas=='true'))+
    geom_point(data=tdf %>% filter(meas!='true'))+facet_wrap(~meas)+
    scale_color_manual(values=cmap_cust)+ theme(legend.position='none')
  plt = plt + labs(x=element_text('time'),
                   y=element_text('signal'))
  plt = plt + scale_x_continuous(breaks=c(0,1))+
    scale_y_continuous(limits = c(-4,4),
                       breaks=c(-2,0,2))
  plt = plt + theme( strip.background=element_blank(),
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
p1=as.ggplot(motivation_fig(1,F,Nm=12,Nm1=NaN))
p2=as.ggplot(motivation_fig(6,T,Nm=12,Nm1=NaN)) 

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
Amps = c(sqrt(2))
Nmvec  = c(12)#c(10,20,30)
acros=seq(0,2*pi,length.out=Nacro)
acros=acros[1:(length(acros)-1)]
freqs  = seq(0,8,length.out=Nfreq)
pars   = expand.grid(Amp=Amps,
                     Nmeas=Nmvec,
                     freq=freqs,
                     acro=acros,
                     type=c('equispaced','optimal'))

pars$power = c(1:dim(pars)[1]) %>% mclapply(mc.cores=8,function(ii){
  x=pars[ii,]
  Amp   = x$Amp
  Nmeas = x$Nmeas
  acro  = x$acro
  freq  = x$freq
  if (x$type=='equispaced'){
    mt = c(1:Nmeas)/Nmeas-1/Nmeas
  }else{
    mt = tirr
  }
  return(evalExactPower(mt,param=list(Amp=Amp,acro=acro,freq=freq)))
})

pars$power = as.numeric(pars$power)

pars =pars %>% mutate(A=Amp,N=Nmeas)

rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),expression(pi/2),
            expression(pi),expression(3*pi/2),
            expression(2*pi))

p3 = pars  %>% ggplot(aes(y=freq,x=acro,fill=power))+
  geom_raster()+
  facet_wrap(~type,ncol=1)+
  scale_x_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  scale_y_continuous(limits=c(0,max(freqs)),
                     breaks =seq(0,max(freqs),4),
                     labels =seq(0,max(freqs),4))+
  scale_fill_viridis_c(limits=c(0,1))+
 # annotate("text", x = 0, y = 6, label = "x", color = "red", size = 8,data=pars[pars$type=='equispaced',]) +
#  geom_abline(slope=0,intercept=6,color='white',linewidth=1.5)+
#  geom_abline(slope=0,intercept=1,color='white',linewidth=1.5)+
  labs(x='acrophase',y='frequency (cycles/day)')
p3

tight_plt=F
if (tight_plt){
  pars %>% filter(type=='optimal') %>% ggplot(aes(y=freq,x=acro,fill=power))+
    geom_raster()+
    facet_wrap(~type,ncol=1)+
    scale_x_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
    scale_y_continuous(limits=c(5.5,6.5),
                       breaks =seq(0,max(freqs),4),
                       labels =seq(0,max(freqs),4))+
    scale_fill_viridis_c(limits=c(0.6,.9))+
    labs(x='acrophase',y='frequency (cycles/day)')
}
plt_width=6

p3 = p3 + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25),
  legend.spacing.x = unit(0, "mm"),
  legend.spacing.y = unit(0, "mm"),
  legend.margin = margin(0, 0, 0, 0)
)
p1=p1&theme(text=element_text(size=fsize))
p2=p2&theme(text=element_text(size=fsize))
p3=p3&theme(text=element_text(size=fsize))



# add raw time plot
Nm=12
tunif   = c(1:Nm)/Nm-1/Nm
tirr    = tau[Xraw[3,]>1e-12]
df=rbind(data.frame(time=tunif,type='equispaced'),
         data.frame(time=tirr,type='optimal'))
dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
head(df)
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point()+
geom_rect(data=dat_bands,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
          inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~type,
             ncol=2)
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,4))
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
p0=plt

Fig = (((p0/p1/p2) + plot_layout(heights=c(.75,6,6)))|p3) +
  plot_annotation(tag_levels='A')+
  plot_layout(widths=c(3.5,1))
show_temp_plt(Fig,6,5)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f1_motivation.png'),
       Fig,
       width=6,height=5,
       device='png',
       dpi=600)
