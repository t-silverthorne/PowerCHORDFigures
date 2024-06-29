#' Evaluate exact power of harmonic regression hypothesis test.
#' 
#' @description
#' Returns exact power of harmonic regression hypothesis test. 
#' 
#' 
#' @param t measurement schedule
#' @param param$Amp amplitude of signal 
#' @param param$freq frequency of signal
#' @param param$acro phase of signal in radians
#' @param alpha type I error, by default \code{alpha=.05}
#' 
#' @note
#' Assumes the noise has mean zero and unit standard deviation.
#' 
#' @return statistical power 
#' @author Turner Silverthorne
#' @export
evalExactPower <- function(t,param,alpha=.05,method='schur',lambda_in=NULL){
# return power of one-frequency cosinor model
  Amp    = param[['Amp']]
  freq   = param[['freq']]
  acro   = param[['acro']]
  N      = length(t)
   
  if (method=='full'){
    A      = matrix(c(0,0,1,0,0,1),nrow=3,byrow=T)
    X      = matrix(c(rep(1,N),cos(2*pi*freq*t),sin(2*pi*freq*t)),ncol=3)
    B      = t(A)%*%solve(t(X)%*%X,A)
    beta   = matrix(c(Amp*cos(acro),Amp*sin(acro)),nrow=2)
    lambda = t(beta)%*%solve(B,beta)
  }else if (method=='schur'){
    A      = matrix(c(0,0,1,0,0,1),nrow=3,byrow=T)
    Xr     = matrix(c(cos(2*pi*freq*t),sin(2*pi*freq*t)),ncol=2)
    D      = t(Xr)%*%Xr
    b      = matrix(c(sum(cos(2*pi*freq*t)),sum(sin(2*pi*freq*t))),ncol=1)
    invB   = D - b%*%t(b)/N
    beta   = matrix(c(Amp*cos(acro),Amp*sin(acro)),nrow=2)
    lambda = t(beta)%*%invB%*%beta
  }else if (method=='ncp'){
    lambda = lambda_in
  }else if (method=='old'){
    cvec   = Amp*cos(2*pi*freq*t-acro)
    lambda = as.numeric(t(cvec)%*%cvec)
  }else{
    stop('unknown method')
  }
  
  f0     = qf(p=1-alpha,df1=2,df2=N-3)
  return(1 - pf(q=f0,df1=2,df2=N-3,ncp=lambda))
}

