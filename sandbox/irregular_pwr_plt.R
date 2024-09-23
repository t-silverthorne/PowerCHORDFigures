mt=tirr
Nacro     = 2^6+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]
fmin=1
fmax=6
pars=expand.grid(acro=acros,
                 freq=seq(5.5,6.5,.01))

pdf=c(1:dim(pars)[1]) %>% lapply(function(ii){
  acro  = pars[ii,]$acro

  param=list(Amp=sqrt(2),freq=pars[ii,]$freq,acro=acro)  
  pwr=evalExactPower(mt,param) 
  return(data.frame(cbind(pars[ii,],data.frame(pwr=pwr))))
}) %>% rbindlist() %>% data.frame()

pdf %>% ggplot(aes(x=acro,y=freq,fill=pwr))+geom_raster()+
  scale_fill_viridis_c()

#pdf %>% filter(freq %in% c(1,6)) %>% 
#  ggplot(aes(x=acro,y=pwr))+geom_line()+facet_wrap(~freq)