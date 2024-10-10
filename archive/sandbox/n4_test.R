require(ggplot2)
require(devtools)
require(data.table)
require(dplyr)
devtools::load_all()
N     = 4
mt    = c(1:N)/N-1/N
acros = seq(0,2*pi,2*pi/2^9)
df=c(1:length(acros)) %>% lapply(function(ii){
  param=list(Amp=10,acro=acros[ii],freq=1)
  data.frame(acro=acros[ii],power=evalExactPower(mt,param,method='full'))
}) %>% rbindlist() %>% data.frame()

df %>% ggplot(aes(x=acro,y=power))+geom_line()