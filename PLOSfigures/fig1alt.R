# setup
source('PLOSfigures/clean_theme.R')
rm()
gc()
n      = 48
tau    = c(1:n)/n -1/n
Xraw   = read.csv2('PLOSfigures/data/cutsdp_sols.csv',header = F,sep=',')
Amp    = sqrt(2) 
Nm     = 12 
tunif  = c(1:Nm)/Nm-1/Nm # equispaced
tirr   = tau[Xraw[3,]>1e-12] # good irregular

Nm1    = 8   
Nm2    = Nm-Nm1
lamb   = 4/24
t1     = c(1:Nm1)/Nm1 -1/Nm1
t2     = c(1:Nm2)/Nm2 -1/Nm2
tirr_b = c(t1*lamb,lamb +t2*(1-lamb))


Nmc     = 1e5
mt      = tunif
acrovec = 2*pi*runif(Nmc)
p24     = .5
pnull   = 0 
p4      = .5
state   = sample(x=c('circ','null','ultra'),
                 prob = c(p24,pnull,p4),
                 size=Nmc,replace=T)

run_sim = function(mtloc,acrovec,state){
  # generate simulated data
  Ydat = matrix(rnorm(Nmc*Nm),nrow=Nmc)
  Ydat[state=='circ',] =Ydat[state=='circ',]+ 
    Amp*cos(outer(acrovec[state=='circ'],2*pi*mtloc,'-'))
  Ydat[state=='ultra',] =Ydat[state=='ultra',]+ 
    Amp*cos(outer(acrovec[state=='ultra'],2*pi*6*mtloc,'-'))
  
  # lomb scargle analysis 
  df = c(1:Nmc) |> mclapply(mc.cores=mc_cores,function(ind){
    lomb = lsp(Ydat[ind,],times=mtloc,plot=F,from = 1,to=8,ofac=10)
    return(data.frame(ind=ind,
                      freq  = lomb$scanned,
                      power = lomb$power))
  }) |> rbindlist() |> data.frame()
  
  # run cosinor at peaks 
  rowc_c = matrixTests::row_cosinor(Ydat,mtloc,1) 
  rowc_u = matrixTests::row_cosinor(Ydat,mtloc,1/6) 
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
resi   = run_sim(tirr,acrovec,state)
resb   = run_sim(tirr_b,acrovec,state)

dfunif = cbind(resu$lombdf,type='equispaced')
dfirr  = cbind(resi$lombdf,type='good alt')
dfirb  = cbind(resb$lombdf,type='bad alt')

dfall  = rbind(dfunif,dfirr,dfirb)

pdf=rbind(cbind(resu$pvdf,type='equispaced'),
          cbind(resi$pvdf,type='good alt'),
          cbind(resb$pvdf,type='bad alt'))

dfgrp = dfall |> group_by(type,freq) |> 
  summarise(mean_power = mean(power))


###########################
# Plot Lomb-Scargle 
###########################
p1 = dfgrp |> ggplot(aes(x=freq,y=mean_power,group=type,color=type))+geom_line()
p1=p1+theme(legend.position='bottom')

###########################
# Plot acro distribution 
###########################
p2a = pdf |>  filter(pval<.05 & type=='equispaced') |> 
  ggplot(aes(x=acro))+geom_histogram(bins=30)+
  facet_wrap(~sig,scales='free_y',nrow=2)
p2a

p2b = pdf |>  filter(pval<.05 & type=='bad alt') |> 
  ggplot(aes(x=acro))+geom_histogram(bins=30)+
  facet_wrap(~sig,scales='free_y',nrow=2)
p2b

p2c = pdf |>  filter(pval<.05 & type=='good alt') |> 
  ggplot(aes(x=acro))+geom_histogram(bins=30)+
  facet_wrap(~sig,scales='free_y',nrow=2)
p2c

p2 = p2a|p2b|p2c


# assemble plot
(p0/p1/p2) + plot_layout(heights=c(1,2,4)) + plot_annotation(tag_levels='A')

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
