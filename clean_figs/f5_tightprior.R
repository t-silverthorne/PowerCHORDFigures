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
tdf = tdf |> filter(window %in% c(12,16,20,24))
dim(tdf)
# naive designs 
ndf = c(12,16,20,24) |> lapply(function(ii){
  N     = 8  
  Nm1   = N-1
  ww    = ii/48
  unm   = c(1:Nm1)/Nm1 - 1/Nm1 
  tvec  = c(ww/2,ww + (1-ww)*unm)
  tvec  = (tvec+0.5) %% 1
  return(data.frame(window=ii,time=tvec,type='naive'))
}) |> rbindlist() |> data.frame()

tdf=rbind(tdf,ndf)


#data.frame(time=2*pi*tvec) |> ggplot(aes(x=time,y=1))+
#  geom_point()+coord_polar(theta='x')+
#  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
#                     labels = c('0','6hr','12hr','18hr'),
#                     limits = c(0,2*pi))



tdf$wlen = 2*pi*tdf$window/48
tdf_full =tdf
tdf = tdf %>%filter(window %in% seq(12,24,4)) 

color1     = "lightblue" 
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
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f3_tightprior1.png'),
       p1,
       width=6,height=2,
       device='png',
       dpi=600)

