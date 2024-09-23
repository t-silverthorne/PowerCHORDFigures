source("figs/fig_settings.R")
devtools::load_all()

if (pub_qual){
  Nmc       = 1e4
  freq_vals = seq(1,30,.05)
}else{
  Nmc       = 1e3
  freq_vals = seq(1,30,.2)
}
mc_cores  = 12 
pars      = expand.grid(freq=freq_vals,
                         Nmeas=c(32,40,48),
                         Amp = c(1),
                         p_osc = c(0.5),
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

plt=df %>% 
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
plt  