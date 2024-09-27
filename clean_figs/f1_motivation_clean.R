source('clean_figs/clean_theme.R')


# load in cutsdp solutions
n    = 48
tau  = c(1:n)/n -1/n
#Xraw = read.csv2('matlab/LMI_formulation/cutsdp_sols.csv',header = F,sep=',')
Xraw = read.csv2('clean_figs/data/cutsdp_sols.csv',header = F,sep=',')
head(Xraw)

Nmc=1e4

Amp     = 2
Nm      = 12 
tunif   = c(1:Nm)/Nm-1/Nm
tirr    = tau[Xraw[3,]>1e-12]
acrovec = 2*pi*runif(Nmc)

###########################
# Acrophase histogram 
###########################
pdf = c(1,6) |> lapply(function(freq){
  Yunif   = Amp*cos(outer(acrovec,2*pi*freq*tunif,'-'))+matrix(rnorm(length(tunif)*Nmc),nrow=Nmc)
  Yirr    = Amp*cos(outer(acrovec,2*pi*freq*tirr,'-'))+matrix(rnorm(length(tirr)*Nmc),nrow=Nmc)

  df_unif = data.frame(meas='equispaced',acro=acrovec,pval= matrixTests::row_cosinor(Yunif,tunif,1/freq) %>% {.$pvalue}) 
  df_irr  = data.frame(meas='optimal',acro=acrovec,pval= matrixTests::row_cosinor(Yirr,tirr,1/freq) %>% {.$pvalue}) 

    
  df_unif = df_unif[df_unif$pval<.05,]
  df_irr  = df_irr[df_irr$pval<.05,]

  df_all  = data.frame(meas='true',acro=acrovec,pval=NA)
  df      = rbind(df_unif,df_irr,df_all)

  df$meas = factor(df$meas,levels=c('true','equispaced','optimal'))
  df$freq=freq
  return(df)
}) |> rbindlist() |> data.frame()

pdf$cmap_var = paste0(pdf$meas,pdf$freq)
pdf = pdf |> mutate(per_label = ifelse(freq==1,'T = 24 hr','T = 4 hr'))
cmap_cust = c('true1'=rgb(.05,0.5,.06),
          'optimal1'=rgb(.05,0.5,.06),
          'equispaced1'=rgb(.05,0.5,.06),
          'true6'=rgb(.05,0.5,.06),
          'optimal6'=rgb(.05,0.5,.06),
          'equispaced6'=rgb(.81,.06,.13))

rad_brk = c(0,pi,2*pi)
rad_lab = c(expression(0),
            expression(pi),
            expression(2*pi))
plt= pdf %>% ggplot(aes(x=acro,fill=cmap_var))+geom_histogram(aes(y=after_stat(density)))+
  facet_grid(per_label~meas)+
  scale_fill_manual(values=cmap_cust)
plt = plt + scale_x_continuous(limits=c(0,2*pi),
                                breaks =rad_brk,
                                labels = rad_lab)
plt = plt+clean_theme()
plt = plt + theme(legend.position='none')
plt = plt + labs(x=element_text('acrophase'),
                  y=element_text('density'))
pbot = plt
p1=plt



###########################
# Raw meas times 
###########################
Nm=12
tunif   = c(1:Nm)/Nm-1/Nm
tirr    = tau[Xraw[3,]>1e-12]
df=rbind(data.frame(time=tunif,type='equispaced'),
         data.frame(time=tirr,type='optimal'))
dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
head(df)
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point()+
geom_rect(data=dat_bands,
  aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
  inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~type,ncol=2)
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,4))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
p0          =plt


###########################
# Heatmap 
###########################
Nfreq = 2^9
Nacro = 2^6+1

Amps = c(sqrt(2))
Nmvec  = c(12)
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
p3 = p3 + clean_theme()
p3 = p3 + theme(
  legend.spacing.x = unit(0, "mm"),
  legend.spacing.y = unit(0, "mm"),
  legend.margin = margin(0, 0, 0, 0)
)

require(patchwork)
Fig = ((p0 / p1) + plot_layout(heights=c(1,6)) | p3)  + 
  plot_layout(widths=c(5,1)) + plot_annotation(tag_levels='A')
Fig
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f1_motivation.png'),
       Fig,
       width=6,height=2.5,
       device='png',
       dpi=600)
