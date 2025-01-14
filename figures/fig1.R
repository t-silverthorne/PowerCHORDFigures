# setup
source('PLOSfigures/clean_theme.R')
rm()
gc()
n      = 48
tau    = c(1:n)/n -1/n
Xraw   = read.csv2('PLOSfigures/data/cutsdp_sols.csv',header = F,sep=',')
Nm     = 12 
tunif  = c(1:Nm)/Nm-1/Nm # equispaced
tirr   = tau[Xraw[3,]>1e-12] # good irregular

Nm1    = 8   
Nm2    = Nm-Nm1
lamb   = 4/24
t1     = c(1:Nm1)/Nm1 -1/Nm1
t2     = c(1:Nm2)/Nm2 -1/Nm2
tirr_b = c(t1*lamb,lamb +t2*(1-lamb))


Nmc     = 1e4
acrovec = 2*pi*runif(Nmc)

Amp     =2 
acrovec = 2*pi*runif(Nmc)
state   = c(rep('circ',Nmc/2),rep('ultra',Nmc/2))

dfplt=c(1:3) |> lapply(function(jj){
  dtype=NaN
  mtloc=NaN
  if (jj==1){
    dtype='equi'
    mtloc=tunif
  }else if (jj==2){
    dtype='irr'
    mtloc=tirr
  }else if (jj==3){
    dtype='irrb'
    mtloc=tirr_b
  }
  Ydat = matrix(rnorm(Nmc*Nm),nrow=Nmc)
  Ydat[state=='circ',] =Ydat[state=='circ',]+ 
    Amp*cos(outer(acrovec[state=='circ'],2*pi*mtloc,'-'))
  Ydat[state=='ultra',] =Ydat[state=='ultra',]+ 
    Amp*cos(outer(acrovec[state=='ultra'],2*pi*6*mtloc,'-'))
  
  # run cosinor for a list of frequencies
  freqs = seq(0.5,7,length.out=1e2)
  
  dfloc = freqs |> mclapply(mc.cores=12,function(freq){
    df_sum = matrixTests::row_cosinor(Ydat,mtloc,1/freq) |> 
      mutate(acro=acrophase*freq*2*pi,
             bin=cut(acro,
                     breaks = seq(0, 2 * pi, length.out = 2^5),
                     include.lowest = TRUE)) |> 
      group_by(bin) |> 
      filter(pvalue<.05) |> 
      summarise(count = n())
    return(cbind(df_sum,data.frame(freq=freq,design=dtype)))    
  }) |> rbindlist() |> data.frame()
}) |> rbindlist() |> data.frame()

#dfplt$bin <- factor(dfplt$bin, levels = unique(dfplt$bin))
dfplt$bin <- factor(dfplt$bin, labels = gsub("\\(|\\[|,.*", "", levels(dfplt$bin)))

dfplt$design =factor(dfplt$design,
                     levels=c('equi','irrb','irr'),
                     labels=c('equispaced','fast-slow','methodically balanced')
                     )

# Create the raster plot
Fig = dfplt |> ggplot( aes(y = freq, x = bin, fill = log10(count)))+
  geom_raster() + facet_wrap(~design,nrow=1)+ 
  scale_fill_viridis_c(name = "Count") +  
  labs( y = "frequency (cycles/day)") +
  scale_x_discrete(
    name = "acrophase (rad)",
    breaks=c(levels(dfplt$bin)[1],
             levels(dfplt$bin)[15],
             levels(dfplt$bin)[31]),
    labels= c("0", "π", "2π")
  )+clean_theme() 
  

ggsave(paste0('PLOSfigures/',
              'fig1spec.png'),
       Fig,
       width=6,height=3.55,
       device='png',
       dpi=600)

# plot acrophase distributions after filtering for significance as a raster plot
