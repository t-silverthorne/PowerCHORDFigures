source('figs/fig_settings.R')
require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(devtools)
require(ggplot2)
require(data.table)
devtools::load_all()

opt_col =rgb(0,.66,.47)
std_col =rgb(.33,.41,.47)

# load in matlab results
ff='matlab/necklace/analyze_sols/optimal_nt_48_.csv'
df = read.csv2(ff,sep=',',header = F)
df = na.omit(df)
bmat = df[,4:dim(df)[2]]
rdf  = df[,1:3]
names(rdf) = c('Nmeas','Nnight','spacing')
am48 = annmatrix(bmat,rann=rdf)

ff='matlab/necklace/analyze_sols/optimal_nt_72_.csv'
df = read.csv2(ff,sep=',',header=F)
df = na.omit(df)
bmat = df[,4:dim(df)[2]]
rdf  = df[,1:3]
names(rdf) = c('Nmeas','Nnight','spacing')
am72 = annmatrix(bmat,rann=rdf)

###############################
# Plot optimal/alt solutions
###############################
nt = 48 
N  = 10 
m  = 2 
if (nt==48){
  am = am48
}else{
  am = am72
}

tau = c(1:nt)/nt - 1/nt
filt = am@''$Nmeas==N & am@''$Nnight==m & am@''$spacing==nt 
x    = as.numeric(am[filt,])
topt = tau[x>1e-12]
topt = (topt-topt[2])%%1

t1   = (1:N)/N - 1/N
t2   = (1:m)/m-1/m
t1   = 0.5*t1
t2   = 0.5+0.5*t2
talt = c(t1,t2)
tdf = rbind(data.frame(t=talt,type='standard'),
          data.frame(t=topt,type='optimal'))
tdf$t = 2*pi*tdf$t
tdf$type = factor(tdf$type,levels=c('standard','optimal'))
plt = tdf %>%  ggplot(aes(x=t,y=1,color=type))+geom_point(show.legend=F)+
  coord_polar(theta='x')+theme_minimal()+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
                     labels = c('0hr','6hr    ','12hr','    18hr'),
                     limits = c(0,2*pi))+
  scale_y_continuous(labels=c())+facet_wrap(~type,nrow=2,strip.position='left')+
  scale_color_manual(values=c('standard'=std_col,'optimal'=opt_col))
plt
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)
plt=plt+theme(text=element_text(size=fsize),legend.position = 'none')
p1=plt
p1

###############################
# Plot power curve 
###############################
rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),
            expression(pi/2),
            expression(pi),
            expression(3*pi/2),
            expression(2*pi))
Nacro     = 2^6+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]
pars=expand.grid(Amp=c(2),
                 Nmeas=N+m,
                 acro=acros,
                 freq=1,
                 type=c('alt design','opt design'))

pwr_vec = c(1:dim(pars)[1]) %>% lapply(function(ii){
  # unpack
  x      = pars[ii,]
  Nmeas  = as.numeric(x[['Nmeas']])
  acro   = as.numeric(x[['acro']])
  freq   = as.numeric(x[['freq']])
  Amp    = as.numeric(x[['Amp']]) 
  param =list(Amp=Amp,freq=freq,acro=acro) 
  
  if(x$type == 'alt design'){
    mt = talt 
  }else if(x$type =='opt design'){
    mt = topt
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
pars$power = as.numeric(pwr_vec)
pars$type = factor(pars$type,levels=c('alt design','opt design'),
                   labels=c('standard','optimal'))
plt=pars %>% ggplot(aes(x=acro,y=power,group=type,color=type))+
  geom_line()+
  ylim(c(0.5,1))+
  scale_color_manual(values=c('standard'=std_col,'optimal'=opt_col))
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),expression(pi/2),
            expression(pi),expression(3*pi/2),
            expression(2*pi))
plt = plt + scale_x_continuous(limits=c(0,2*pi),
                               breaks =rad_brk,
                               labels = rad_lab)
plt=plt+theme(text=element_text(size=fsize),legend.position='bottom')
p2=plt


###############################
# Plot power recovery 
###############################
gen_df = function(am){
  c(1:dim(am)[1]) %>% lapply(function(ii){
    nt   = dim(am)[2]
    tau = c(1:nt)/nt - 1/nt
    x    = as.numeric(am[ii,])
    topt = tau[x>1e-12]
    
    N    = am@''[ii,]$Nmeas
    m    = am@''[ii,]$Nnight
    t1   = (1:N)/N - 1/N
    t2   = (1:m)/m-1/m
    t1   = 0.5*t1
    t2   = 0.5+0.5*t2
    talt = c(t1,t2)

    Ntot  = N+m  
    tunif = c(1:Ntot)/Ntot-1/Ntot
    
      
    min_eig_opt  = evalMinEig(topt,1)
    min_eig_alt  = evalMinEig(talt,1)
    min_eig_unif = evalMinEig(tunif,1)
    
    cbind(am@''[ii,],data.frame(ncp_opt=min_eig_opt,
                                ncp_alt=min_eig_alt,
                                ncp_unif=min_eig_unif))
  }) %>% rbindlist() %>% data.frame()
}
df=gen_df(am48)

head(df)
opt_col =rgb(0,.66,.47)
std_col =rgb(.33,.41,.47)


df_opt = df
df_opt$ncp_rel = df_opt$ncp_opt/df_opt$ncp_unif
df_opt$Ntot    = df_opt$Nmeas + df_opt$Nnight
df_opt = df_opt[,c('Ntot','ncp_rel')]
df_opt$type ='optimal'


df_alt= df
df_alt$ncp_rel = df_alt$ncp_alt/df_alt$ncp_unif
df_alt$Ntot    = df_alt$Nmeas + df_alt$Nnight
df_alt = df_alt[,c('Ntot','ncp_rel')]
df_alt$type ='standard'

df_tall = rbind(df_opt,df_alt)

plt=df_tall %>% ggplot(aes(x=Ntot,y=ncp_rel,color=type,fill=type))+
  geom_point(show.legend=F)+geom_line(show.legend=F)+
  scale_color_manual(values=c('standard'=std_col,'optimal'=opt_col))
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+theme(text=element_text(size=fsize))
p3=plt

Fig=(p1+p2+p3)+ plot_annotation(tag_levels='A')+
  plot_layout(widths=c(.5,1,1),guides='collect') & theme(legend.position='bottom') 

show_temp_plt(Fig,6,2.75)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f3_tightprior.png'),
       Fig,
       width=6,height=2.75,
       device='png',
       dpi=600)