source('figures/fig_settings.R')
###########################
# bad scaling
###########################
Nm   = 8
if (pub_qual){
  Nmc  = 1e4
  nrep = 5 
}else{
  Nmc  = 1e3
  nrep = 3 
}
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
cmap_cust = c('MC'=rgb(0,.55,.55),
              'exact'='blue',
              'approximate'='orange')
plt=pdf %>% ggplot(aes(x=scale,y=Power,group=method,color=method))+geom_line()+geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper))+
  scale_color_manual(values=cmap_cust)

plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)

# might need
plt = plt + labs(x=element_text('scale'),
                 y=element_text('power'))


# always need
plt=plt+theme(text=element_text(size=fsize))

p1=plt

###########################
# Bias in wrong formula
###########################
freq  = 1
Nm    = 24 
if (pub_qual){
  Nmc   = 1e4
  nrep  = 1e4
}else{
  Nmc   = 1e3
  nrep  = 1e3
}

param = list(Amp=1,freq=1,acro=0)

set.seed(1)
df = c(1:nrep) %>% mclapply(mc.cores=8,function(ii){
  tvec = runif(Nm)
  pwr_mc      = evalMonteCarloPower(tvec,param,Nmc)
  pwr_exact   = evalExactPower(tvec,param)
  pwr_approx  = evalExactPower(tvec,param,method='old')
  return(data.frame(pwr_mc=pwr_mc,pwr_exact=pwr_exact,pwr_approx=pwr_approx))
}) %>% rbindlist() %>% data.frame()

df$bias_exact  = df$pwr_exact-df$pwr_mc
df$bias_approx = df$pwr_approx -df$pwr_mc

pdf = data.frame(bias=c(df$bias_exact,df$bias_approx),method=c(rep('exact',nrep),rep('approx',nrep)))

pdf$method = factor(pdf$method,levels=c('exact','approx'),labels=c('exact','approximate'))
cmap_cust=c('exact'='blue',
            'approximate'='orange')
plt=pdf %>% ggplot(aes(x=bias,fill=method,group=method))+
  geom_histogram(position='identity',bins=100)+facet_wrap(~method,nrow=1)+
  scale_fill_manual(values=cmap_cust)
plt

plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)

plt = plt + labs(x=element_text('MC power difference'))
plt=plt+theme(text=element_text(size=fsize))

p2=plt
Fig = p2/p1+plot_annotation(tag_levels='A')

show_temp_plt(Fig,6,4)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'fig2.png'),
       Fig,
       width=6,height=4,
       device='png',
       dpi=600)

#df1 = df[,c('pwr_mc','pwr_approx')]
#names(df1)=c("pwr_mc",'pwr')
#df2 = df[,c('pwr_mc','pwr_exact')]
#names(df2)=c("pwr_mc",'pwr')
#
#df1$method = 'approx'
#df2$method = 'exact'
#
#df=rbind(df1,df2)
#
#p2=df %>% ggplot(aes(x=pwr_mc,y=pwr,group=method,color=method))+geom_point(size=.5)+geom_abline(slope=1,intercept=0)
