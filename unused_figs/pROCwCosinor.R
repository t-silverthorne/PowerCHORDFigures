# Check amplitude has expected effect
require(pROC)
require(data.table)
require(dplyr)
require(parallel)
require(ggplot2)
p_osc = 0.5
Nm    = 12
Nmc   = 1e5
mt    = c(1:Nm)/Nm - 1/Nm

freq  = 1

check_AUC = function(Amp){
  Ydat  = matrix(rnorm(Nmc*Nm),nrow=Nmc)
  state = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+
        Amp*cos(outer(2*pi*runif(N_osc),2*pi*freq*mt,'-'))
  
  pvec = matrixTests::row_cosinor(Ydat,mt,period=1/freq) |> (\(x) x$pvalue)()
  
  return(pROC::roc(as.numeric(state=='osc'),pvec,direction='>') |> (\(x) x$auc )() )
}
check_AUC(1)
check_AUC(2)


# Check intuition with cosinor model
fvec =seq(1,5,.5)
Amp=1.5
check_cosinor = function(ftrue){
  Ydat  = matrix(rnorm(Nmc*Nm),nrow=Nmc)
  state = sample(c('osc','non_osc'),Nmc,replace = T,c(p_osc,1-p_osc))
  N_osc = sum(state=='osc')
  Ydat[state=='osc',] = Ydat[state=='osc',]+
        Amp*cos(outer(2*pi*runif(N_osc),2*pi*ftrue*mt,'-'))
  
  pvdf = fvec |> lapply(function(freq){
    pvec = matrixTests::row_cosinor(Ydat,mt,period=1/freq) |> (\(x) x$pvalue)()
    return(data.frame(
      idx=c(1:dim(Ydat)[1]),
      pvalue = pvec,
      freq = freq
      ))
  }) |> rbindlist() |> data.frame() |> 
    group_by(idx) |> 
    summarize(pvalue=min(pvalue))
  pvec = pvdf$pvalue
  return(pROC::roc(as.numeric(state=='osc'),pvec,direction='>') |> (\(x) x$auc )() )
}

true_freqs = seq(1,5,.05)
adf = true_freqs |> mclapply(mc.cores=10,function(ftrue){
  auc = check_cosinor(ftrue) |> as.numeric() 
  return(data.frame(auc=auc,ftrue=ftrue))
  }) |> rbindlist() |> data.frame()

adf |> ggplot(aes(x=ftrue,y=auc))+geom_line()