#'Use Monte Carlo method to estimate power 
#' 
#' @description
#' Performs harmonic regression on a simulated dataset of independent samples
#' and returns the portion of samples that have statistically significant p-values.
#' Useful for comparison with the exact expression for statistical power,
#' see [eval_exact_power]. Relies on [rowCosinor].
#' 
#' @param tvec vector of measurement times
#' @param param$Amp amplitude of signal 
#' @param param$freq frequency fo signal
#' @param param$acro phase of signal in radians
#' @param Nmc number of Monte Carlo samples
#' @param alpha type I error, default value \code{alpha=.05} 
#' 
#' @return Monte Carlo estimate of power
#' @author Turner Silverthorne
#' @export
evalMonteCarloPower<-function(tvec,param,Nmc,alpha=.05,method='Ftest',Nperm=1e2){
  Amp  = param$Amp
  freq = param$freq
  acro = param$acro
  Ydat = replicate(Nmc,{Amp*cos(2*pi*freq*tvec -acro) + rnorm(length(tvec))}) %>% t
  if (method=='Ftest'){
    pwr  = rowCosinor(Ydat,tvec,1/freq) %>% {.$pvalue <alpha} %>% mean()
  }else if(method=='perm'){
    #stop('perm not implemented')
    pvec_true = rowCosinor(Ydat,tvec,1/freq) %>% {.$pvalue}
    perm_mat  = matrix(rep(NaN,Nmc*Nperm),nrow=Nmc,ncol = Nperm)
    for (perm in c(1:Nperm)){
      Pdat = matrix(rep(NaN,Nmc*length(tvec)),nrow=Nmc,ncol = length(tvec))
      for (ii in c(1:Nmc)){
        Pdat[ii,] = Ydat[ii,sample(1:length(tvec))]
      }
      perm_mat[,perm]=rowCosinor(Pdat,tvec,1/freq) %>% {.$pvalue}
    } 
    pmat_true = replicate(Nperm,pvec_true)
    perm_mat < pmat_true
    R = (perm_mat<pmat_true) %>% rowSums()
    pwr =mean((R+1)/(Nperm+1) < alpha)
  }else{
    stop('Unknown method')
  }
  return(pwr)
}
