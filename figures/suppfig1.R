source('PLOSfigures/clean_theme.R')
###########################
# bad scaling
###########################
Nm   = 8   # sample size
Nmc  = 1e4 # number of monte carlo samples
nrep = 5   # number of replicates for error bars

mt   =  c(1:Nm)/Nm -1/Nm
Amp  = 2
freq = 1
acro = pi

# compute Monte Carlo Power
sc= seq(.1,1,.05)
scales = rep(sc,nrep)
df=c(1:length(scales)) %>% lapply(function(ii){
  return(data.frame(scale=scales[ii],
                    power=evalMonteCarloPower(mt*scales[ii],
                                              Amp=Amp,freq=freq,acro=acro,Nmc=Nmc)))
}) %>% rbindlist() %>% data.frame()


head(df)
df_grp = df %>% group_by(scale) %>% summarize(Power=mean(power),
                                     lower=quantile(power)[2],
                                     upper=quantile(power)[4])

df_grp$method = 'MC'

# compute power using general closed-form expression
df_exact = c(1:length(sc)) %>% lapply(function(ii){
  return(data.frame(scale=scales[ii],
                    Power=evalPower(mt*scales[ii],Amp=Amp,freq=freq,acro=acro)))
}) %>% rbindlist() %>% data.frame()

df_exact$lower = NaN
df_exact$upper = NaN
df_exact$method='general'

# compute power using closed-form expression which only holds for equispaced
df_wrong = c(1:length(sc)) %>% lapply(function(ii){
  return(data.frame(scale=scales[ii],
                    Power=evalPower(mt*scales[ii],Amp=Amp,freq=freq,acro=acro,method='equispaced')))
}) %>% rbindlist() %>% data.frame()

df_wrong$lower = NaN
df_wrong$upper = NaN
df_wrong$method= 'equispaced'

pdf=rbind(df_grp,df_exact,df_wrong)

pdf$method = factor(pdf$method,levels=c('MC','general','equispaced'))
cmap_cust = c('MC'=rgb(0,.55,.55),
              'general'='blue',
              'equispaced'='orange')

# compare the three power estimates in plot
plt=pdf %>% ggplot(aes(x=scale,y=Power,group=method,color=method))+
  geom_line()+geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper))+
  scale_color_manual(values=cmap_cust)
plt=plt+clean_theme()
plt = plt + labs(x=element_text('scale'),
                 y=element_text('power'))
plt = plt + theme(legend.position='bottom')
p1=plt

###########################
# Bias in wrong formula
###########################
freq  = 1
Nm    = 24 
Nmc   = 1e4
nrep  = 1e4
Amp   = 1
freq  = 1
acro  = 0

set.seed(1)
df = c(1:nrep) %>% mclapply(mc.cores=mc_cores,function(ii){
  tvec = runif(Nm)
  pwr_mc      = evalMonteCarloPower(tvec,Amp=Amp,freq=freq,acro=acro,Nmc)
  pwr_exact   = evalPower(tvec,Amp=Amp,freq=freq,acro=acro)
  pwr_approx  = evalPower(tvec,
                               Amp=Amp,freq=freq,acro=acro,
                               method='equispaced')
  return(data.frame(pwr_mc=pwr_mc,pwr_exact=pwr_exact,pwr_approx=pwr_approx))
}) %>% rbindlist() %>% data.frame()

df$bias_exact  = df$pwr_exact-df$pwr_mc
df$bias_approx = df$pwr_approx -df$pwr_mc

pdf = data.frame(bias=c(df$bias_exact,df$bias_approx),method=c(rep('general',nrep),rep('equispaced',nrep)))

pdf$method = factor(pdf$method,levels=c('general','equispaced'),labels=c('general','equispaced'))
cmap_cust=c('general'='blue',
            'equispaced'='orange')
plt=pdf %>% ggplot(aes(x=bias,fill=method,group=method))+
  geom_histogram(position='identity',bins=100)+facet_wrap(~method,nrow=2)+
  scale_fill_manual(values=cmap_cust)
plt=plt+clean_theme()
plt = plt + theme(legend.position='bottom')
plt = plt + labs(x=element_text('MC power difference'))
p2=plt
p2
Fig =(p2|p1)+plot_annotation(tag_levels='A') + plot_layout(widths=c(1,1))

ggsave('PLOSfigures/suppfig1.png',
       Fig,
       width=6,height=3,
       device='png',
       dpi=600)