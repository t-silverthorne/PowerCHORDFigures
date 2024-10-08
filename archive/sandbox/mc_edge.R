require(devtools)
require(data.table)
require(ggplot2)
devtools::load_all()

# measurement schedule
N = 8
mt    = c(1:N)/N-1/N
mt    = c(mt,mt)
param = list(Amp=1,freq=1,acro=pi/2)

# compare exact and monte carlo
paste0('Exact power   : ',evalExactPower(mt,param))
paste0('Monte Carlo   : ', evalMonteCarloPower(mt,param,1e5))

# get mininum eigenvalue
cat('noncentrality: ', evalMinEig(mt,1))


# check N=3 edge case
N     = 3
mt    = c(1:N)/N-1/N
acros = seq(0,2*pi,2*pi/2^6)
df=c(1:length(acros)) %>% lapply(function(acro){
  param=list(Amp=3,acro=acro,freq=1)
  data.frame(acro=acro,power=evalExactPower(mt,param,method='full'))
}) %>% rbindlist() %>% data.frame()

df %>% ggplot(aes(x=acro,y=power))+geom_line()
