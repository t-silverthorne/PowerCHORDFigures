source("clean_figs/clean_theme.R")
Nmeas  = 32 

acro_dist = 'average'

Nmc    = 1e4
freq   = 1 
Amp    = 2 
p_osc  = 0.5 

sols   = readRDS('clean_figs/data/powerCHORD_even_sols.RDS')
filt = sols@Nmeas==Nmeas & sols@fmin == 1 & sols@fmax==Nmeas/2 & sols@method=='diffEVCR'
mt_alt = sols[filt,]
mt_alt = as.numeric(mt_alt)
mt_alt = mt_alt[!is.nan(mt_alt)]
mt_alt = sort(mt_alt)

mt_unif = c(1:Nmeas)/Nmeas - 1/Nmeas

Ydat   = matrix(rnorm(Nmc*Nmeas),nrow=Nmc)
state  = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
N_osc  = sum(state=='osc')

for (mt in list(mt_unif,mt_alt)){
  if (acro_dist =='average'){
    Ydat[state=='osc',] = Ydat[state=='osc',]+
      Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*mt,'-'))
  }else if (acro_dist =='worst'){
    acros     = seq(0,2*pi,length.out=2^8+1)
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
  print(roc$auc)
}