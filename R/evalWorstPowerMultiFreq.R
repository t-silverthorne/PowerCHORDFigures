evalWorstPowerMultiFreq=function(mt,param,alpha=.05,method='eig',returnType='min'){
  Amp   = param[['Amp']]
  fmin  = param[['fmin']]
  fmax  = param[['fmax']]
  Nfreq = param[['Nfreq']]

  freqs = seq(from=fmin,to=fmax,length.out=Nfreq)

  worst_pwrs=freqs %>% sapply(function(freq){
    ploc=list(Amp=Amp,freq=freq)
    evalWorstPower(mt,ploc,alpha,method)
  })

  if (returnType=='min'){
    return(min(unlist(worst_pwrs)))
  }else if(returnType=='all'){
    return(unlist(worst_pwrs))
  }else{
    stop('Unknown return type')
  }
}