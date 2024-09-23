source('figs/fig_settings.R')
require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(devtools)
require(ggplot2)
require(data.table)
devtools::load_all()

# load matlab solutions
ff='matlab/necklace/night_filt2/window_sols.csv'
df = read.csv2(ff,sep=',',header = F)
df = na.omit(df)

wvec=seq(12,24,2)
n = 48
tau = c(1:n)/n-1/n

###########################
# polar plot of raw solutions 
###########################
wstart = rep(NaN,dim(df)[1])

for (ind in c(1:dim(df)[1])){
  bvec = df[ind,] %>% as.matrix() %>% as.numeric()
  lenw = wvec[ind]
  pv   = rep(1,lenw)
  found_window = F
  for (ii in c(0:(n-1))){
    if (t(bvec[1+(0:(lenw-1) + ii)%%n])%*%pv == 1){
      if (!found_window){
        found_window=T
        wstart[ind]= ii
      }
    }
  }
}
wstart

# rotate dataframe so that windows align
for (ii in c(1:dim(df)[1])){
  rvec = df[ii,c(1 + ( 0:(n-1) + wstart[ii] -1  + n/2)%%n )]
  df[ii,] = rvec
}

# dataframe of measurement times
tdf = c(1:dim(df)[1]) %>% lapply(function(ii){data.frame(window=wvec[ii],time=tau[df[ii,]>0])}) %>% 
  rbindlist() %>% data.frame()

tdf$wlen = 2*pi*tdf$window/48
tdf_full =tdf
tdf = tdf %>% filter(window!=18)


color1     = "lightblue" 
color2     = "darkred"  
color_func = colorRamp(c(color1, color2))
num_colors = 6 # Number of colors to generate
f_colors   = rgb(color_func(seq(0, 1, length.out = num_colors)), maxColorValue = 255)

setNames(f_colors, unique(tdf$wlen))

tdf$window = factor(tdf$window,unique(tdf$window),paste0(unique(tdf$window/2),'hr rest'))
plt = tdf %>% filter(window!=18) %>%  ggplot(aes(x=2*pi*time,y=1))+
  geom_point()+facet_wrap(~window,nrow=3,dir='v')+
  geom_col(data=tdf,aes(y = 0.1, x = pi,fill=as.factor(wlen)),just = 0,width=tdf$wlen) +  
  coord_polar(theta='x')+theme_minimal()+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
                     labels = c('0','6hr','12hr','18hr'),
                     limits = c(0,2*pi))+
  scale_y_continuous(labels=c())+
  scale_fill_manual(values = f_colors)
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  axis.text.x = element_text(vjust = 0.25),
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)
plt=plt+theme(text=element_text(size=fsize),legend.position = 'none')
p1=plt


# acrophase plot
tdf=tdf_full
N=8
tunif = c(1:N)/N - 1/N
tcstr1 = tdf %>% filter(window==12) %>% select(time) %>% as.matrix() %>% as.numeric() 
tcstr2 = tdf %>% filter(window==24) %>% select(time) %>% as.matrix() %>% as.numeric() 

Nacro     = 2^6+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]

pars=expand.grid(acro=acros,type=c('equispaced','constrained1','constrained2'))

pdf = c(1:dim(pars)[1]) %>% lapply(function(ii){
  x=pars[ii,]
  acro = as.numeric(x[['acro']]) 
  freq = 1
  Amp  = 2 
  param=list(Amp=Amp,freq=freq,acro=acro)
  if (x[['type']]=='equispaced'){
    mt=tunif 
  }else if(x[['type']]=='constrained1'){
    mt=tcstr1
  }else if(x[['type']]=='constrained2'){
    mt=tcstr2
  }else{
    stop('unknown type')
  }
  return(data.frame(cbind(x,power=evalExactPower(mt,param))))
}) %>% rbindlist() %>%  data.frame()

rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),
            expression(pi/2),
            expression(pi),
            expression(3*pi/2),
            expression(2*pi))
cmap_manual = c('equispaced'='black','constrained1'='lightblue',constrained2='darkred')
lmap_manual = c('equi'='dashed','constr'='solid') 
pdf$is_constrained='equi'
pdf[pdf$type%in%c('constrained1','constrained2'),]$is_constrained='constr'
plt=pdf %>% ggplot(aes(x=acro,y=power,group=type,color=type,linetype=is_constrained))+
  geom_line()+scale_linetype_manual(values=lmap_manual)+
  scale_color_manual(values=cmap_manual)+
  scale_x_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])
plt=plt+theme(legend.position='none')
plt=plt+labs(x=element_text('acrophase (rad)'),
                 y=element_text('power'))
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),text=element_text(size=fsize)
)
p2=plt

# spread plot
tdf=tdf_full
wvec 
rel_ncp=c(1:length(wvec)) %>% lapply(function(ii){
rel_ncp=tdf %>% filter(window==wvec[ii]) %>% select(time) %>% as.matrix() %>% as.numeric() %>% 
  evalMinEig(freq=1)/4
data.frame(rel_ncp=as.numeric(rel_ncp))
}) %>% rbindlist()

pdf=cbind(data.frame(window=wvec),rel_ncp)
plt=pdf %>% ggplot(aes(x=window/2,y=rel_ncp))+geom_line()+geom_point()+
  labs(x='rest duration (hr)',y='relative noncentrality\n parameter')+
  scale_x_continuous(breaks=seq(12,24,4)/2)
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),text=element_text(size=fsize)
)
p3=plt


Fig = (p1 |(p2/p3) + plot_layout(heights=c(1,1)) ) + plot_annotation(tag_levels='A')

show_temp_plt(Fig,6,4)


ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f3_tightprior.png'),
       Fig,
       width=6,height=4,
       device='png',
       dpi=600)