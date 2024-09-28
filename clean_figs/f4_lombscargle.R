source("clean_figs/clean_theme.R")
pub_qual=T
if (pub_qual){
  Nmc       = 1e3 
  freq_vals = seq(1,30,.05)
}else{
  Nmc       = 10 
  freq_vals = seq(1,30,10)
}

mc_cores  = 8 
pars      = expand.grid(freq=freq_vals,
                         Nmeas=c(32,40,48),
                         Amp = c(1),
                         p_osc = c(0.5),
                         type=c('irregular','equispaced'))

#sols   = readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
sols   = readRDS('clean_figs/data/powerCHORD_even_sols.RDS')
dim(pars)


df=c(1:dim(pars)[1]) %>% mclapply(mc.cores=mc_cores,function(ind){
  freq        = pars[ind,]$freq
  Amp         = pars[ind,]$Amp
  p_osc       = pars[ind,]$p_osc
  Nmeas       = pars[ind,]$Nmeas
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
  Ydat[state=='osc',] = Ydat[state=='osc',]+Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*mt,'-'))
  
  # simualte p-values
  pvdf = c(1:dim(Ydat)[1]) %>% lapply(function(ii){
    x              = Ydat[ii,]
    lomb_std       = lsp(x,times=mt,plot=F,normalize = 'standard')
    return(data.frame(p_method='std',pval =lomb_std$p.value,state=state[ii]))
  }) %>% rbindlist() %>% data.frame()
  
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
saveRDS(df,'clean_figs/data/lomb_out.RDS')
#plt=df %>% 
#  ggplot(aes(x=freq,y=AUC,group=type,color=type))+geom_line()+
#  geom_vline(aes(xintercept = Nmeas / 2), linetype = "dashed", color = "black")+
#  facet_grid(Nmeas~Amp)+theme(legend.position='bottom')+labs(x='frequency (cycles/day)')
#plt = plt+clean_theme()
#plt=plt+theme(legend.position='bottom',
#              legend.direction = "horizontal")
#plt  
