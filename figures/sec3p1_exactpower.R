fsize=9
require(dplyr)
require(data.table)
require(ggplot2)
require(ggplotify)
require(patchwork)
require(devtools)
load_all()
theme_set(theme_classic()) 

Nm   = 8
Nmc  = 1e4
nrep = 5
mt = c(1:Nm)/Nm -1/Nm
param = list(Amp=2,freq=1,acro=pi)

sc= seq(.1,1,.05)
scales = rep(sc,nrep)
df=c(1:length(scales)) %>% lapply(function(ii){
  return(data.frame(scale=scales[ii],
                    power=evalMonteCarloPower(mt*scales[ii],param,Nmc=Nmc)))
}) %>% rbindlist() %>% data.frame()


head(df)
df_grp = df %>% group_by(scale) %>% summarize(Power=mean(power),
                                     lower=quantile(power)[2],
                                     upper=quantile(power)[4])

df_grp$method = 'MC'


df_exact = c(1:length(sc)) %>% lapply(function(ii){
  return(data.frame(scale=scales[ii],
                    Power=evalExactPower(mt*scales[ii],param)))
}) %>% rbindlist() %>% data.frame()

df_exact$lower = NaN
df_exact$upper = NaN
df_exact$method='exact'

df_wrong = c(1:length(sc)) %>% lapply(function(ii){
  return(data.frame(scale=scales[ii],
                    Power=evalExactPower(mt*scales[ii],param,method='old')))
}) %>% rbindlist() %>% data.frame()

df_wrong$lower = NaN
df_wrong$upper = NaN
df_wrong$method= 'approximate'

pdf=rbind(df_grp,df_exact,df_wrong)

pdf$method = factor(pdf$method,levels=c('MC','exact','approximate'))
p1=pdf %>% ggplot(aes(x=scale,y=Power,group=method,color=method))+geom_line()+geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper))



freq  = 1
Nm    = 24 
Nmc   = 1e3
nrep  = 1e3
param = list(Amp=1,freq=1,acro=0)
d1    = replicate(nrep,{runif(Nm) %>% {evalExactPower(.,param)}})
d2    = replicate(nrep,{runif(Nm) %>% {evalExactPower(.,param,method='old')}})
d3    = replicate(nrep,{runif(Nm) %>% {evalMonteCarloPower(.,param,Nmc=Nmc)}})

d1  = data.frame(method='exact',power=d1)
d2  = data.frame(method='approximate',power=d2)
d3  = data.frame(method='MC',power=d3)
pdf = rbind(d1,d2,d3)

pdf$method = factor(pdf$method,levels=c('MC','exact','approximate'))
pdf %>% ggplot(aes(x=power,group=method,fill=method)) + geom_histogram( position="identity")+
  facet_wrap(~method,ncol=1)




