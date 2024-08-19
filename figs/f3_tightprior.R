source('figs/fig_settings.R')
require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(devtools)
require(ggplot2)
require(data.table)
devtools::load_all()

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
nt = 72 
N  = 6
m  = 2 
if (nt==48){
  am = am48
}else{
  am = am72
}

tau = c(1:nt)/nt - 1/nt
filt = am@''$Nmeas==N & am@''$Nnight==m & am@''$spacing==nt 
x    = as.numeric(am[filt,])
topt = tau[x>0]

t1   = (1:N)/N - 1/N
t2   = (1:m)/m-1/m
t1   = 0.5*t1
t2   = 0.5+0.5*t2
talt = c(t1,t2)
data.frame(t=talt) %>% ggplot(aes(x=t,y=1))+geom_point()+xlim(c(0,1))+
  geom_vline(xintercept = .5)
  x

plt = data.frame(t=talt*2*pi) %>% ggplot(aes(x=t,y=1))+geom_point()+
  coord_polar(theta='x')+theme_minimal()+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
                     labels = c('0','6hr','12hr','18hr'),
                     limits = c(0,2*pi))+
  scale_y_continuous(labels=c())
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+theme(text=element_text(size=fsize))
p1=plt

plt = data.frame(t=((topt-topt[2])%%1)*2*pi) %>% ggplot(aes(x=t,y=1))+geom_point()+
  coord_polar(theta='x')+theme_minimal()+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
                     labels = c('0','6hr','12hr','18hr'),
                     limits = c(0,2*pi))+
  scale_y_continuous(labels=c())
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+theme(text=element_text(size=fsize))
p2=plt

###########################
# power heatmap 
###########################
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
plt=pars %>% ggplot(aes(x=acro,y=power,group=type,color=type))+geom_line()+ylim(c(0,1))
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+theme(text=element_text(size=fsize))
p3=plt


gen_df = function(am){
  c(1:dim(am)[1]) %>% lapply(function(ii){
    nt   = dim(am)[2]
    tau = c(1:nt)/nt - 1/nt
    x    = as.numeric(am[ii,])
    topt = tau[x>0]
    
    N    = am@''[ii,]$Nmeas
    m    = am@''[ii,]$Nnight
    t1   = (1:N)/N - 1/N
    t2   = (1:m)/m-1/m
    t1   = 0.5*t1
    t2   = 0.5+0.5*t2
    talt = c(t1,t2)
    
    min_eig_opt = evalMinEig(topt,1)
    min_eig_alt = evalMinEig(talt,1)
    
    cbind(am@''[ii,],data.frame(ncp_opt=min_eig_opt,
                                ncp_alt=min_eig_alt))
  }) %>% rbindlist() %>% data.frame()
}
df = rbind(gen_df(am48),gen_df(am72))

df$dncp = 100*(df$ncp_opt - df$ncp_alt)/df$ncp_alt

plt=df %>% ggplot(aes(x=Nmeas+2,y=dncp,group=spacing,shape=as.factor(spacing)))+geom_point()+
  geom_line()+labs(x='total measurement budget',y='ncp improvement (%)',color='grid')
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)

# always need
plt=plt+theme(text=element_text(size=fsize))

p4 = plt

((p1|p2)/p3/p4) + plot_annotation(tag_levels='A')