source('figures/clean_theme.R')
Nfreq    = 2^10
Amp      = 1

# load differential evolution results 
am = readRDS('figures/data/diffEvolveOutput.RDS')
am = am[am@''$method=='diffEVCR',]
df = am@''

###########################
# Plot of raw solutions 
###########################
fmax       = 24
fmin       = 1
Nmeas_vals = c(24,32,48)

sloc     = am[am@''$Nmeas %in% Nmeas_vals & am@''$fmin==fmin & am@''$fmax==fmax,]
sloc[1,] = sloc[1,]-min(sloc[1,],na.rm=T)
sloc[2,] = sloc[2,]-min(sloc[2,],na.rm=T)
splt     = sloc %>% stack()
splt     = splt[!is.nan(splt$value),]
splt     = splt %>% mutate(N=Nmeas)

plt=splt[splt$Nmeas %in% Nmeas_vals, ] %>% ggplot(aes(x=value,y=0))+geom_point(size=.5)+
  facet_wrap(~N,nrow=3,strip.position='right')
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(labels=seq(0,24,4),
  breaks=seq(0,1,4/24),
  limits=c(0,1))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
psol=plt


###########################
# Alt plot of raw solutions 
###########################
set.seed(2)
mtopt = am[am@''$Nmeas ==24 & am@''$fmin==fmin & am@''$fmax==fmax,]
mtopt = mtopt %>% as.numeric()
mtopt = mtopt[!is.nan(mtopt)]
mtopt = mtopt-min(mtopt)
sc = 15/60/24

mtunif   = c(1:24)/24 -1/24
mtopt_j  = mtopt + rnorm(24,0,sc)
mtunif_j = mtunif + rnorm(24,0,sc)

tdf = rbind(data.frame(time=mtunif,type='equispaced',jitter='raw'),
  data.frame(time=mtunif_j,type='equispaced',jitter='jittered'),
  data.frame(time=mtopt,type='irregular',jitter='raw'),
  data.frame(time=mtopt_j,type='irregular',jitter='jittered')
  )
tdf$jitter=factor(tdf$jitter,c('raw','jittered'))
head(tdf)

plt = tdf %>%  ggplot(aes(x=time,y=0))+geom_point(size=.5)+
  facet_grid(jitter~type)
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(labels=seq(0,24,4),
  breaks=seq(0,1,4/24),
  limits=c(0,1))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
psol=plt
psol

###########################
# Robustness to measurement timing 
###########################
Nmvec  = c(24,32,48)
nrep   = 100
scales = c(0,seq(1,30,2))/60/24
pars   = expand.grid(scale=scales,type=c('irregular','equispaced'),Nm=Nmvec)
sols   = readRDS('figures/data/diffEvolveOutput.RDS')


df=c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ii){
  sc   = pars[ii,]$scale
  type = pars[ii,]$type
  Nm   = pars[ii,]$Nm
  fmax_loc = Nm/2
  fmin = 1
  Amp = 1
  fmax=fmax_loc
  Nfreq=Nfreq
  if (type=='equispaced'){
    mt = c(1:Nm)/Nm -1/Nm 
  }else{
    filt = sols@Nmeas==Nm & sols@fmin == fmin & sols@fmax==fmax_loc & sols@method=='diffEVCR'
    mt = sols[filt,]
    mt = as.numeric(mt)
    mt = mt[!is.nan(mt)]
    if(length(mt)!=Nm){
      stop('wrong length meas vec')
    }
  }
  cbind(pars[ii,],
        data.frame(power=replicate(nrep,{
          evalWorstPowerMultiFreq(mt+rnorm(Nm,0,sd=sc),
                                  fmin=fmin,
                                  fmax=fmax,
                                  Nfreq=Nfreq,
                                  Amp=Amp)})))
}
) %>% rbindlist() %>% data.frame()

df_grp = df %>% group_by(scale,type,Nm) %>% summarize(Power=mean(power),
                                                      lower=quantile(power)[2],
                                                      upper=quantile(power)[4])
df_grp$time = df_grp$scale*24*60

df_grp

# summary statistic: first scale where you cross halfway point
thresh_df = df_grp  |>
  group_by(Nm,scale) |> 
  summarise(dpower = Power[type=='irregular'] - 
              Power[type=='equispaced' ],.groups='drop')
tdf0=thresh_df |> filter(scale==0) |> select(Nm,dpower) |> rename(dpower0=dpower)

thresh_df_summary = thresh_df |> left_join(tdf0,by="Nm") |> 
  filter(dpower < dpower0/2) |> 
  group_by(Nm) |> 
  slice_min(scale) |> 
  mutate(time=scale*24*60)

thresh_df_summary

# summary statistic: delta between noiseless and final
tds_2 = df_grp |> filter(type=='irregular') |> 
  group_by(Nm) |> 
  summarise(dpower = Power[scale==0] - Power[scale==max(scale)],
            max_scale = max(scale))
tds_2

plt=df_grp  |>  filter(scale>0) |>  ggplot(aes(x=time,y=Power,group=type,color=type)) +geom_line() +
  geom_errorbar(aes(ymin=lower,ymax=upper))+
  facet_wrap(~Nm,nrow=1)
plt = plt+clean_theme()
plt = plt + labs(x=element_text('collection time standard deviation (minutes)'),
                 y=element_text('power'),
                 color='design')
plt = plt + theme(legend.position='right')
prob = plt

Fig=psol/prob + plot_layout(heights=c(2,1.5))+plot_annotation(tag_levels='A')

ggsave(
  filename = "vector_figures/Main05.pdf",
  plot = Fig,
  device = "pdf",
  width = 6,
  height = 2.75,
  units = "in"    
)