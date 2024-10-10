source('clean_figs/clean_theme.R')
# load in cutsdp solutions
N     = 12
n     = 48
tau   = c(1:n)/n -1/n
freqs = c(2,4,6,8,10,12)
#Xraw=read.csv2('matlab/LMI_formulation/cutsdp_sols.csv',header = F,sep=',')
Xraw=read.csv2('clean_figs/data/cutsdp_sols.csv',header = F,sep=',')
head(Xraw)

###########################
# Plot of raw solutions
###########################
topt  = tau[as.numeric(Xraw[6,])>0]
tunif = c(1:N)/N - 1/N
length(topt)==length(tunif)

df = rbind(data.frame(time=topt,type='optimal'),
      data.frame(time=tunif,type='equispaced'))

dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point()+
  geom_rect(data=dat_bands,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
            inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~type,
             ncol=1,strip.position='top')
plt=plt+clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,2))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
p1=plt
###########################
# ncp across freqs 
###########################
fvec = seq(0,16,.01)
pars = expand.grid(freq=fvec,type=c('optimal','equispaced'))

edf = c(1:dim(pars)[1]) %>% lapply(function(ii){
  x=pars[ii,]
  freq=x[['freq']]
  if(x[['type']]=='equispaced'){
    mt = tunif
  }else if(x[['type']]=='optimal'){
    mt = topt 
  }else{
    stop('unknown type')
  }
  return(cbind(pars[ii,],data.frame(emin =evalMinEig(mt,freq))))
}) %>% rbindlist() %>% data.frame()

color_scale= c('optimal'=rgb(0.13, 0.67, 0.8),
               'equispaced'=rgb(.43,.21,.1))
plt=edf %>% ggplot(aes(x=freq,y=emin,group=type,color=type))+geom_line()+
  geom_vline(xintercept = 12,linetype='dashed')+
  geom_vline(xintercept = 1,linetype='dashed')+
  scale_x_continuous(limits=c(0,16),breaks=seq(0,16,4))+
  scale_color_manual(values=color_scale)+
  labs(x='frequency (cycles/day)',y='noncentrality parameter')
plt = plt+clean_theme()
plt = plt+theme(legend.position='bottom')
p2=plt

###########################
# power across acro 
###########################
#Nacro = 2^6+1
#acros = seq(0,2*pi,length.out=Nacro)
#acros = acros[1:(length(acros)-1)]
#pars = expand.grid(acro=acros,freq=c(1,12))
#pdf = c(1:dim(pars)[1]) %>% lapply(function(ii){
#  x=pars[ii,]
#  freq=as.numeric(x[['freq']])
#  acro=as.numeric(x[['acro']])
#  param=list(Amp=2,acro=acro,freq=freq)
#  return(cbind(pars[ii,],data.frame(power=evalExactPower(topt,param))))
#}) %>% rbindlist() %>% data.frame()
#pdf %>% ggplot(aes(x=acro,y=power))+geom_line()+facet_wrap(~freq)

data.frame(time = (topt %% (1/12))) %>% ggplot(aes(x=time,y=1))+geom_point()


###########################
# Pareto frontier  
###########################

nrep = 1e4

#random designs
rdf=c(1:nrep) %>% lapply(function(ii){
  mt = runif(N)  
  return(data.frame(type='random',
                    eig1=evalMinEig(mt,1),
                    eig12=evalMinEig(mt,12)))
}) %>% rbindlist() %>% data.frame()
# equispaced
rdf=rbind(rdf,data.frame(type='equispaced',eig1=evalMinEig(tunif,1),eig12=evalMinEig(tunif,12)),
data.frame(type='optimal',eig1=evalMinEig(topt,1),eig12=evalMinEig(topt,12)))

alpha_scale = c('optimal'=1,'equispaced'=1,'random'=.1)
size_scale  = c('optimal'=2.5,'equispaced'=2.5,'random'=.7)
type_scale = c('optimal'=18,'equispaced'=18,'random'=19)
color_scale= c('optimal'=rgb(0.13, 0.67, 0.8),
               'equispaced'=rgb(.43,.21,.1),
               'random'='black')
plt=rdf %>% ggplot(aes(x=eig1,y=eig12,color=type,size=type,
                       alpha=type,shape=type))+
  geom_point()+
  scale_size_manual(values=size_scale)+
  scale_alpha_manual(values = alpha_scale)+
  scale_color_manual(values=color_scale)+
  scale_shape_manual(values=type_scale)+
  labs(x='noncentrality (f=1)',y='noncentrality (f=12)')
plt=plt+clean_theme()
plt = plt+theme(legend.position='bottom')
p3=plt
p3

Fig= (p1/( p2+p3 )) + plot_layout(heights=c(1,3))+plot_annotation(tag_levels='A')
show_temp_plt(Fig,6,3.5)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f_cutsdp1.png'),
       Fig,
       width=6,height=3.5,
       device='png',
       dpi=600)
# optimal

