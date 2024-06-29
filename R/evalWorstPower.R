evalWorstPower=function(mt,param,alpha=.05,method='eig'){
  Amp    = param[['Amp']]
  freq   = param[['freq']]
  N      = length(t)
  if (length(freq)>1){
    stop('Use evalWorstPowerMutliFreq for handling multiple frequencies')
  }
  if (method=='eig'){
    A       = matrix(c(0,0,1,0,0,1),nrow=3,byrow=T)
    Xr      = matrix(c(cos(2*pi*freq*mt),sin(2*pi*freq*mt)),ncol=2)
    D       = t(Xr)%*%Xr
    b       = matrix(c(sum(cos(2*pi*freq*mt)),sum(sin(2*pi*freq*mt))),ncol=1)
    invB    = D - b%*%t(b)/length(mt)
    ncp     = eigen(invB) %>% {.$values} %>% min() %>% {.*Amp^2} 
    min_pwr = evalExactPower(mt,param,method='ncp',lambda_in=ncp)
    min_pwr 
  }else if(method=='test'){
    Nacro = 2^12
    acro = seq(0,2*pi,length.out=Nacro+1)
    acro = acro[1:Nacro]
    min_pwr= acro %>% sapply(function(phi){
      ploc      = param
      ploc$acro = phi
      return(evalExactPower(mt,ploc,alpha))
    }) %>% min()
    min_pwr
  }else{
    stop('unknown method')
  }
  return(min_pwr)
}