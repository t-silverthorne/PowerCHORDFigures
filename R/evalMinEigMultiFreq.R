evalMinEigMultiFreq=function(mt,param,alpha=.05){
  fmin  = param[['fmin']]
  fmax  = param[['fmax']]
  Nfreq = param[['Nfreq']]

  freqs = seq(from=fmin,to=fmax,length.out=Nfreq)

  min_eigs=freqs %>% sapply(function(freq){
    evalMinEig(mt,freq)
  })
  min(unlist(min_eigs))
}