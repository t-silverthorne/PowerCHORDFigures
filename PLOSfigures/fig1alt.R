# setup
source('plosfigures/clean_theme.r')
rm()
gc()
n      = 48
tau    = c(1:n)/n -1/n
xraw   = read.csv2('plosfigures/data/cutsdp_sols.csv',header = f,sep=',')
amp    = sqrt(2) 
nm     = 12 
tunif  = c(1:nm)/nm-1/nm # equispaced
tirr   = tau[xraw[3,]>1e-12] # good irregular

nm1    = 8   
nm2    = nm-nm1
lamb   = 4/24
t1     = c(1:nm1)/nm1 -1/nm1
t2     = c(1:nm2)/nm2 -1/nm2
tirr_b = c(t1*lamb,lamb +t2*(1-lamb))


nmc     = 1e5
mt      = tunif
acrovec = 2*pi*runif(nmc)
p24     = .5
pnull   = 0 
p4      = .5
state   = sample(x=c('circ','null','ultra'),
                 prob = c(p24,pnull,p4),
                 size=nmc,replace=t)

run_sim = function(mtloc,acrovec,state){
  # generate simulated data
  ydat = matrix(rnorm(nmc*nm),nrow=nmc)
  ydat[state=='circ',] =ydat[state=='circ',]+ 
    amp*cos(outer(acrovec[state=='circ'],2*pi*mtloc,'-'))
  ydat[state=='ultra',] =ydat[state=='ultra',]+ 
    amp*cos(outer(acrovec[state=='ultra'],2*pi*6*mtloc,'-'))
  
  # lomb scargle analysis 
  df = c(1:nmc) |> mclapply(mc.cores=mc_cores,function(ind){
    lomb = lsp(ydat[ind,],times=mtloc,plot=f,from = 1,to=8,ofac=30)
    return(data.frame(ind=ind,
                      freq  = lomb$scanned,
                      power = lomb$power))
  }) |> rbindlist() |> data.frame()
  
  # run cosinor at peaks 
  rowc_c = matrixtests::row_cosinor(ydat,mtloc,1) 
  rowc_u = matrixtests::row_cosinor(ydat,mtloc,1/6) 
  pvalc =  rowc_c |> (\(x) x$pvalue)()
  acroc =  rowc_c |> (\(x) 2*pi*x$acrophase*1)()
  pvalu =  rowc_u |> (\(x) x$pvalue)()
  acrou =  rowc_u |> (\(x) 2*pi*x$acrophase*6)()
 

  dfp = rbind(data.frame(acro=acrou,
                   pval=pvalu,
                   state=state,
                   sig ='ultradian',
                   acro_true=acrovec,
                   rfreq=6),
              data.frame(acro=acroc,
                   pval=pvalc,
                   state=state,
                   sig ='circadian',
                   acro_true=acrovec,
                   rfreq=1
                   ))
  return(list(lombdf=df,pvdf=dfp))
}
resu   = run_sim(tunif,acrovec,state)
resb   = run_sim(tirr_b,acrovec,state)

dfunif = cbind(resu$lombdf,type='equispaced')
dfirb  = cbind(resb$lombdf,type='bad alt')

dfall  = rbind(dfunif,dfirb)

pdf=rbind(cbind(resu$pvdf,type='equispaced'),
          cbind(resb$pvdf,type='bad alt'))

dfgrp = dfall |> group_by(type,freq) |> 
  summarise(mean_power = mean(power))

###########################
# plot raw measurement times 
###########################
nm=12
tunif   = c(1:Nm)/Nm-1/Nm
tirr    = tau[Xraw[3,]>1e-12]
df=rbind(data.frame(time=tunif,type='equispaced'),
         data.frame(time=tirr_b,type='fast-slow alternative'))
dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
head(df)
df$type = factor(df$type,levels=c('equispaced','fast-slow alternative','methodically balanced'))
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point(size=.5)+
geom_rect(data=dat_bands,
  aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
  inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~type,ncol=1)
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,4))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
p0  = plt



###########################
# Plot Lomb-Scargle 
###########################

dfgrp$type <- factor(dfgrp$type, 
                     levels = c("equispaced", "bad alt", "good alt"), 
                     labels = c("equispaced", "fast-slow", "methodically balanced"))
custom_colors = c(
  c("equispaced"='black',
    "fast-slow"=rgb(1,0,0.5),
    "methodically balanced"=rgb(.36,.54,.66)
      )
)
plt = dfgrp |> ggplot(aes(x=freq,y=mean_power,group=type,color=type))+
  geom_line()+
  geom_vline(xintercept = 1, linetype = "dashed")+
  geom_vline(xintercept = 6, linetype = "dashed")
  
plt = plt + clean_theme()
plt = plt+theme(legend.position='bottom')
plt = plt + labs(x='frequency (cycles/day)',y='average intensity')
plt = plt + scale_color_manual(values=custom_colors)
plt
p1=plt
p1=p1+theme(legend.position='right')
###########################
# Plot acro distribution 
###########################
p2a = pdf |>  filter(pval<.05 & type=='equispaced') |> 
  ggplot(aes(x=acro_true))+geom_histogram(bins=30)+
  facet_wrap(~sig,scales='free_y',nrow=2)
p2a

p2b = pdf |>  filter(pval<.05 & type=='bad alt') |> 
  ggplot(aes(x=acro_true))+geom_histogram(bins=30)+
  facet_wrap(~sig,scales='free_y',nrow=2)
p2b


p2 = p2a|p2b

###########################
# Nicer plot of acro distribution 
###########################
pdf$typef <- factor(pdf$type, 
                     levels = c("equispaced", "bad alt"), 
                     labels = c("equispaced", "fast-slow alternative"))

pdf$cmap_var = paste0(pdf$typef,pdf$rfreq)
pdf = pdf |> mutate(per_label = ifelse(rfreq==1,'T = 24 hr','T = 4 hr'))
cmap_cust = c(
          'equispaced1'=rgb(.05,0.5,.06),
          'fast-slow alternative1'=rgb(.81,.06,.13),
          'fast-slow alternative6'=rgb(.05,0.5,.06),
          'equispaced6'=rgb(.81,.06,.13))

plt = pdf |> filter(pval<.05) |> 
  ggplot(aes(x=acro_true,fill=cmap_var))+
  geom_histogram()+facet_grid(per_label~typef)+
  scale_fill_manual(values=cmap_cust)
plt = plt + scale_x_continuous(limits=c(0,2*pi),
                                breaks =c(rad_brk[1],rad_brk[5]),
                                labels =c(rad_lab[1],rad_lab[5]))
plt = plt+clean_theme()
plt = plt + theme(legend.position='none')
plt = plt + labs(x=element_text('acrophase'),
                  y=element_text('count'))
plt
###########################
# assemble plot
###########################
p1=p1+theme(legend.position='bottom')
Fig = (((p0|p1)+
          plot_layout(widths = c(1,2))&plot_layout(guides='collect')&
          theme(legend.position='bottom'))/p2) +
  plot_layout(heights=c(1,2.5)) +plot_annotation(tag_levels='A')

show_temp_plt(Fig,6,3.5)
#Fig=(p0/p1/p2) + plot_layout(heights=c(0.5,2,4)) + plot_annotation(tag_levels='A')
ggsave(paste0('PLOSfigures/',
              'fig1.png'),
       Fig,
       width=6,height=3.5,
       device='png',
       dpi=600)

#pdf |> filter(pval<.05 & type=='bad alt' & sig=='circadian' & state=='circ' ) |> 
#    ggplot(aes(x=acro_true))+geom_histogram(bins=30)
#mt = tunif
#fnow=5.95k
#Ydat = matrix(rnorm(1e3*length(mt)),ncol=length(mt))
#r1= row_cosinor(Ydat,mt,period=1/fnow)
#r2= rowCosinor(Ydat,mt,per=1/fnow,method='ginv')
#
#data.frame(qr_phase=r1$acrophase,
#         ginv_phase=r2$phase) |> 
#  ggplot(aes(x=qr_phase,y=ginv_phase))+geom_point(size=.5)
