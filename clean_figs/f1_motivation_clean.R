source('clean_figs/clean_theme.R')


# load in cutsdp solutions
n    = 48
tau  = c(1:n)/n -1/n
#Xraw = read.csv2('matlab/LMI_formulation/cutsdp_sols.csv',header = F,sep=',')
Xraw = read.csv2('clean_figs/data/cutsdp_sols.csv',header = F,sep=',')
head(Xraw)

Nmc=1e5

Amp     = sqrt(2) 
Nm      = 12 
tunif   = c(1:Nm)/Nm-1/Nm # equispaced
tirr    = tau[Xraw[3,]>1e-12] # good irregular

# bad irregular
Nm1    = 8   
Nm2    = Nm-Nm1
lamb   = 4/24
t1     = c(1:Nm1)/Nm1 -1/Nm1
t2     = c(1:Nm2)/Nm2 -1/Nm2
tirr_b = c(t1*lamb,lamb +t2*(1-lamb))


acrovec = 2*pi*runif(Nmc)

###########################
# Acrophase histogram 
###########################
pdf = c(1,6) |> lapply(function(freq){
  Yunif   = Amp*cos(outer(acrovec,2*pi*freq*tunif,'-'))+matrix(rnorm(length(tunif)*Nmc),nrow=Nmc)
  Yirr    = Amp*cos(outer(acrovec,2*pi*freq*tirr,'-'))+matrix(rnorm(length(tirr)*Nmc),nrow=Nmc)
  Yirr_b  = Amp*cos(outer(acrovec,2*pi*freq*tirr_b,'-'))+matrix(rnorm(length(tirr_b)*Nmc),nrow=Nmc)

  df_unif   = data.frame(meas='equispaced',acro=acrovec,pval= matrixTests::row_cosinor(Yunif,tunif,1/freq) %>% {.$pvalue}) 
  df_irr    = data.frame(meas='good alternative',acro=acrovec,pval= matrixTests::row_cosinor(Yirr,tirr,1/freq) %>% {.$pvalue}) 
  df_irr_b  = data.frame(meas='bad alternative',acro=acrovec,pval= matrixTests::row_cosinor(Yirr_b,tirr_b,1/freq) %>% {.$pvalue}) 

    
  df_unif   = df_unif[df_unif$pval<.05,]
  df_irr    = df_irr[df_irr$pval<.05,]
  df_irr_b  = df_irr_b[df_irr_b$pval<.05,]

  #df_all  = data.frame(meas='true',acro=acrovec,pval=NA)
  #df      = rbind(df_unif,df_irr,df_all)
  df      = rbind(df_unif,df_irr,df_irr_b)

  df$meas = factor(df$meas,levels=c('equispaced','bad alternative','good alternative'))
  df$freq=freq
  return(df)
}) |> rbindlist() |> data.frame()

pdf$cmap_var = paste0(pdf$meas,pdf$freq)


pdf = pdf |> mutate(per_label = ifelse(freq==1,'T = 24 hr','T = 4 hr'))
cmap_cust = c('good alternative1'=rgb(.05,0.5,.06),
          'equispaced1'=rgb(.05,0.5,.06),
          'bad alternative1'=rgb(.81,.06,.13),
          'good alternative6'=rgb(.05,0.5,.06),
          'bad alternative6'=rgb(.05,0.5,.06),
          'equispaced6'=rgb(.81,.06,.13))

rad_brk = c(0,pi,2*pi)
rad_lab = c(expression(0),
            expression(pi),
            expression(2*pi))
plt= pdf %>% filter(meas!='true') %>% 
ggplot(aes(x=acro,fill=cmap_var))+
  geom_histogram()+
  facet_grid(per_label~meas)+
  scale_fill_manual(values=cmap_cust)
plt = plt + scale_x_continuous(limits=c(0,2*pi),
                                breaks =rad_brk,
                                labels = rad_lab)
plt = plt+clean_theme()
plt = plt + theme(legend.position='none')
plt = plt + labs(x=element_text('acrophase'),
                  y=element_text('count'))
pbot = plt
p1=plt



###########################
# Raw meas times 
###########################
Nm=12
tunif   = c(1:Nm)/Nm-1/Nm
tirr    = tau[Xraw[3,]>1e-12]
df=rbind(data.frame(time=tunif,type='equispaced'),
         data.frame(time=tirr,type='good alternative'),
         data.frame(time=tirr_b,type='bad alternative'))
dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
head(df)
df$type = factor(df$type,levels=c('equispaced','bad alternative','good alternative'))
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point(size=.5)+
geom_rect(data=dat_bands,
  aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
  inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~type,ncol=3)
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,4))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
p0          =plt


require(patchwork)
Fig = ((p0 / p1) + plot_layout(heights=c(1,7)) )  + 
  plot_layout(widths=c(3,1)) + plot_annotation(tag_levels='A')
show_temp_plt(Fig,6,2.5)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f1_motivation.png'),
       Fig,
       width=6,height=2.5,
       device='png',
       dpi=600)
