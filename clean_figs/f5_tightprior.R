source("clean_figs/clean_theme.R")

# load matlab solutions
#ff='matlab/necklace/night_filt2/window_sols.csv'
ff='clean_figs/data/window_sols.csv'
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
tdf = c(1:dim(df)[1]) %>% lapply(function(ii){
  data.frame(window=wvec[ii],time=tau[df[ii,]>1e-12],type='optimal')}) |> 
  rbindlist() |> data.frame()

# naive designs 
ndf = seq(12,24,2) |> lapply(function(ii){
  N     = 8  
  Nm1   = N-1
  ww    = ii/48
  unm   = c(1:Nm1)/Nm1 - 1/Nm1 
  tvec  = c(ww/2,ww + (1-ww)*unm)
  tvec  = (tvec+0.5) %% 1
  return(data.frame(window=ii,time=tvec,type='naive'))
}) |> rbindlist() |> data.frame()

tdf=rbind(tdf,ndf)

tdf$wlen = 2*pi*tdf$window/48
tdf_full =tdf
tdf = tdf %>%filter(window %in% seq(12,24,4)) 

color1     = rgb(0.36,0.54,.66)
color2     = "darkred"  
color_func = colorRamp(c(color1, color2))
num_colors = 6 # Number of colors to generate
f_colors   = rgb(color_func(seq(0, 1, length.out = num_colors)), maxColorValue = 255)

setNames(f_colors, unique(tdf$wlen))

tdf$window = factor(tdf$window,unique(tdf$window),paste0(unique(tdf$window/2),'hr rest'))
plt = tdf %>% filter(window!=18) %>% 
  ggplot(aes(x=2*pi*time,y=1))+
  geom_point()+facet_grid(type~window)+
  geom_col(data=tdf,aes(y = 0.1, x = pi,fill=as.factor(wlen)),just = 0,width=tdf$wlen) +  
  coord_polar(theta='x',clip='off')+theme_minimal()+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
                     labels = c('0','6','12','18'),
                     limits = c(0,2*pi))+
  scale_y_continuous(labels=c(),limits = c(0,1.2))+
  scale_fill_manual(values = f_colors)+
  labs(x='time (hr)')
plt=plt+clean_theme()
plt = plt + theme(
  plot.margin = margin(0,0,0,0),
  #axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line = element_blank()
)
plt=plt+theme(legend.position = 'none')
p1=plt
p1


###########################
# power as function of acro
###########################
tdf=tdf_full
N=8
tunif = c(1:N)/N - 1/N
tcstr1 = tdf %>% filter(window==12) %>% select(time) %>% as.matrix() %>% as.numeric() 
tcstr2 = tdf %>% filter(window==24) %>% select(time) %>% as.matrix() %>% as.numeric() 

Nacro     = 2^6+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]

pars=expand.grid(acro=acros,window=c(12,24),
                 type=c('optimal','naive'))

pdf = c(1:dim(pars)[1]) %>% lapply(function(ii){
  x=pars[ii,]
  acro = as.numeric(x[['acro']]) 
  freq = 1
  Amp  = 2.5 
  param=list(Amp=Amp,freq=freq,acro=acro)
  mt = tdf %>% filter(window==x[['window']] & type == x[['type']]) %>% 
    select(time) %>% unlist %>% as.numeric()
  return(data.frame(cbind(x,power=evalExactPower(mt,param))))
}) %>% rbindlist() %>%  data.frame()
head(pdf)

cmap_manual = c('12'=rgb(0.36,0.54,.66),'24'='darkred')
lmap_manual = c('naive'='dashed','optimal'='solid') 
pdf$gvar = paste0(pdf$type,pdf$window)
pdf$type = factor(pdf$type,levels=c('optimal','naive'))
plt= pdf %>% ggplot(aes(x=acro,y=power,group=gvar,color=as.factor(window),linetype=type))+
  geom_line()+scale_linetype_manual(values=lmap_manual)+
  scale_color_manual(values=cmap_manual)+
  scale_x_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  guides(color='none',linetype=guide_legend(title=NULL))
plt=plt+clean_theme()
plt=plt+theme(legend.position='bottom',
              legend.direction = "horizontal")
plt=plt+labs(x=element_text('acrophase (rad)'),
                 y=element_text('power'))
p2=plt
p2

p1 / p2


###########################
# relative ncp across windows
###########################
tdf=tdf_full
wvec 

pars = expand.grid(window=wvec,type=c('naive','optimal'))

pdf = c(1:dim(pars)[1]) %>% lapply(function(ii){
  x  = pars[ii,]
  mt = tdf %>% filter(window==x[['window']] & type == x[['type']]) %>% 
    select(time) %>% unlist %>% as.numeric()
  rel_ncp = evalMinEig(t=mt,freq=1)/4
  return(cbind(pars[ii,],data.frame(rel_ncp=rel_ncp)))
}) %>% rbindlist() %>% data.frame()

pdf$type = factor(pdf$type,levels=c('optimal','naive'))
plt=pdf %>% ggplot(aes(x=window/2,y=rel_ncp,group=type,linetype=type))+geom_line()+geom_point()+
  labs(x='rest duration (hr)',y='relative noncentrality\n parameter')+scale_linetype_manual(values=lmap_manual)+
  scale_x_continuous(breaks=seq(12,24,4)/2)+
  guides(linetype=guide_legend(title=NULL))
plt=plt+clean_theme()
plt=plt+theme(legend.position='bottom',
              legend.direction = "horizontal")
p3=plt
p3

Fig = (p1 | (p2/p3)  )  + plot_layout(guides='collect',widths=c(2,1)) +
  plot_annotation(tag_levels='A')& theme(legend.position='bottom') & guides(color='none',fill='none')
show_temp_plt(Fig,6,3)


ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f3_tightprior1.png'),
       Fig,
       width=6,height=3,
       device='png',
       dpi=600)

