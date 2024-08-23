source('figs/fig_settings.R')
if (pub_qual){
  Nacro     = 2^6
  Nmc       = 1e4
  freq_vals = seq(1,30,.05)
}else{
  Nacro     = 2^5
  Nmc       = 5e3
  freq_vals = seq(1,30,.2)
}
mc_cores  = 12 
pars      = expand.grid(freq=freq_vals,
                         Nmeas=c(32,40,48),
                         Amp = c(1.5,2),
                         p_osc = c(0.5),
                         stat_method=c('lomb'),
                         acro_dist = c('worst'),
                         type=c('WCP','equispaced'))

sols   = readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
dim(pars)
df=c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ind){#parallel inside
  freq        = pars[ind,]$freq
  Amp         = pars[ind,]$Amp
  p_osc       = pars[ind,]$p_osc
  Nmeas       = pars[ind,]$Nmeas
  stat_method = pars[ind,]$stat_method
  acro_dist   = pars[ind,]$acro_dist
  type        = pars[ind,]$type
  
  if (type=='equispaced'){
    mt = c(1:Nmeas)/Nmeas -1/Nmeas 
  }else{
    filt = sols@Nmeas==Nmeas & sols@fmin == 1 & sols@fmax==Nmeas/2 & sols@method=='diffEVCR'
    mt = sols[filt,]
    mt = as.numeric(mt)
    mt = mt[!is.nan(mt)]
    if(length(mt)!=Nmeas){
      stop('wrong length meas vec')
    }
  }
  
  # simulate data
  Nmeas               = length(mt)
  Ydat                = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
  state               = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc               = sum(state=='osc')
  if (acro_dist =='average'){
    Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*mt,'-'))
  }else if (acro_dist =='worst'){
    acros     = seq(0,2*pi,length.out=Nacro)
    acros     = acros[1:(length(acros)-1)]
    powers    = acros %>% sapply(function(acro){
      param=list(freq=freq,acro=acro,Amp=Amp)
      evalExactPower(mt,param)
    }) 
    worst_ind = which.min(powers)
    acro = acros[worst_ind]
    
    signal = Amp*cos(2*pi*freq*mt-acro)
    smat   = t(matrix(rep(signal,N_osc),ncol=N_osc))
    Ydat[state=='osc',] = Ydat[state=='osc',]+smat
  }else{
    stop('unknown acro distribution')
  }
  
  # simualte p-values
  if (stat_method=='cosinor'){
    pvdf = data.frame(pval=rowCosinor(Ydat,mt,per=1/freq) %>% {.$pvalue},
                      state=state)
  }else if(stat_method=='lomb'){
    pvdf = c(1:dim(Ydat)[1]) %>% lapply(function(ii){
      x              = Ydat[ii,]
      lomb_std       = lsp(x,times=mt,plot=F,normalize = 'standard')
      return(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]))
    }) %>% rbindlist() %>% data.frame()
  }else{
    stop('unknown stat method')
  }
  
  # FDR correction
  pval   = pvdf$pval
  ostate = pvdf$state
  
  # record AUC
  roc=pROC::roc(as.numeric(ostate=='osc'),pval,direction='>')
  # record FPR and TPR at alpha=.05
  num_P   = sum(ostate=='osc') 
  num_N   = sum(ostate=='non_osc') 
  num_TP  = sum(ostate=='osc'     & pval < .05) 
  num_FP  = sum(ostate=='non_osc' & pval < .05)
  TPR     = num_TP/num_P
  FPR     = num_FP/num_N
  
  return(data.frame(cbind(pars[ind,],data.frame(AUC=roc$auc,TPR=TPR,FPR=FPR))))
}) %>% rbindlist() %>% data.frame()

# ticks
# y label
# borders
# font size

plt=df %>% filter(acro_dist=='worst' ) %>% 
  ggplot(aes(x=freq,y=AUC,group=type,color=type))+geom_line()+
  geom_vline(aes(xintercept = Nmeas / 2), linetype = "dashed", color = "black")+
  facet_grid(Nmeas~Amp)+theme(legend.position='bottom')+labs(x='frequency (cycles/day)')
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+theme(text=element_text(size=9))
p1=plt
p1
show_temp_plt(p1,6,3)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f5_lombscargl.png'),
       p1,
       width=6,height=3,
       device='png',
       dpi=600)

saveRDS(df,'temp_wide.RDS')