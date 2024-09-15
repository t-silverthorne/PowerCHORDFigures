require(devtools)
require(devtools)
Xraw=read.csv2('matlab/LMI_formulation/cutsdp_sols.csv',header = F,sep=',')
head(Xraw)
n=48
tau = c(1:n)/n -1/n
freqs=c(2,4,6,8,10,12)
ii=1

sanity_check=F
if(sanity_check){
  ###########################
  # Verify optimality 
  ###########################
  df=c(1:length(freqs)) %>% lapply(function(ii){
    mt = tau[as.numeric(Xraw[ii,])>0]
    if(length(mt)==12){
      eig1=evalMinEig(mt,freq=1)
      eig2=evalMinEig(mt,freq=freqs[ii])
      data.frame(freq1=1,freq2=freqs[ii],eig1=eig1,eig2=eig2)
    }else{
      stop('wrong length')
    }
  }) %>% rbindlist() %>% data.frame()
  
  
  ###########################
  # Check power is flat 
  ###########################
  Nacro     = 2^6+1
  acros     = seq(0,2*pi,length.out=Nacro)
  acros     = acros[1:(length(acros)-1)]
  pars=expand.grid(acro=acros,
                   index=c(1:6),
                   type=c('fmin','fmax'))
  
  pdf=c(1:dim(pars)[1])%>% lapply(function(ii){
    acro  = pars[ii,]$acro
    index = pars[ii,]$index
    type  = pars[ii,]$type
  
    mt = tau[as.numeric(Xraw[index,])>0]
    if(type=='fmin'){
      param=list(Amp=sqrt(2),freq=1,acro=acro)  
    } else{
      param=list(Amp=sqrt(2),freq=freqs[index],acro=acro)  
    } 
    pwr=evalExactPower(mt,param) 
    return(data.frame(cbind(pars[ii,],data.frame(pwr=pwr))))
  }) %>% rbindlist() %>% data.frame()
    
}

###########################
# Plot of raw solutions
###########################
df=c(1:length(freqs)) %>% lapply(function(ii){
  mt = tau[as.numeric(Xraw[ii,])>0]
  mt = mt-min(mt) 
  data.frame(time=mt,fmax=freqs[ii])
}) %>% rbindlist() %>% data.frame()

dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
df$fmax_lab = paste0('[1,',df$fmax,']')
df$fmax_lab = factor(df$fmax_lab,levels=c('[1,2]','[1,4]','[1,6]',
                                   '[1,8]','[1,10]','[1,12]'))
head(df)
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point()+
geom_rect(data=dat_bands,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
          inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~fmax_lab,
             ncol=1,strip.position='left')
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,2))
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
p2=plt



###########################
# Comparison with trivial bound 
###########################
require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(patchwork)
am=readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
df = am@''

dfgrp1 = df[,c('Nmeas','fmin','fmax','method','ncp')] %>%
  pivot_wider(names_from = method, values_from = ncp, names_prefix = "ncp_") 
dfgrp2 = df[,c('Nmeas','fmin','fmax','method','upper')] %>%
  pivot_wider(names_from = method, values_from = upper, names_prefix = "upper_") 
dfgrp=merge(dfgrp1,dfgrp2)


plt=dfgrp %>% ggplot(aes(x=ncp_diffEV,y=Nmeas/2,color=fmax,shape=as.factor(Nmeas)))+geom_point()+
  geom_abline(slope=1,intercept=0)+
  scale_x_continuous(trans='log2')+ 
  scale_y_continuous(trans='log2') +scale_color_viridis_c(option='plasma')
plt = plt+scale_y_continuous(breaks=seq(8,50,8))+
  scale_x_continuous(breaks=seq(4,24,4))
plt = plt + labs(x='differential evolution score',
                 y='Nmeas/2',
                 color='maximum frequency',
                 shape='budget')
plt=plt+guides(color=guide_colorbar(title.position='bottom',
                                    position='bottom'))
plt=plt+guides(shape=guide_legend(title.position='bottom',
                                  position = 'bottom'))
plt
plt=plt+theme(text=element_text(size=fsize),
              legend.title = element_text(hjust=0.5)
)
plt
p1=plt
Fig1 = p1
show_temp_plt(Fig1,6,2)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f_cutsdp1.png'),
       Fig1,
       width=6,height=2,
       device='png',
       dpi=600)


Fig2 = p2 
show_temp_plt(Fig2,6,3)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f_cutsdp2.png'),
       Fig2,
       width=6,height=3,
       device='png',
       dpi=600)