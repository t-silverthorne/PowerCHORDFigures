# check how equispaced vs random designs perform for signals that violate model assumptions
require(devtools)
require(ggplot2)
require(data.table)
source('figures/clean_theme.R')
require(PowerCHORD)

n    = 12
tt   = (0:n)/n
tt   = tt[1:n]

Nsamp = 1e4

Amin = 1
Amax = 3
racro = function(){runif(1,0,2*pi)}
rAmp  = function(){runif(1,Amin,Amax)}
sim_noise             = function(tt){rnorm(length(tt))}
sim_osc_noise         = function(tt){rnorm(length(tt))*(1+rAmp()*cos(2*pi*tt - racro()))}
sim_cosinor           = function(tt,noise=T){rAmp()*cos(2*pi*tt - racro())+noise*rnorm(length(tt))}
sim_cosinor_osc_amp   = function(tt,noise=T){(  (1+runif(1,0,.95)*cos(2*pi*tt-racro()))*rAmp()  )*cos(2*pi*tt - racro())+noise*rnorm(length(tt))}
sim_cosinor_osc_noise = function(tt,noise=T){rAmp()*cos(2*pi*tt - racro())+noise*rnorm(length(tt))*(1+rAmp()*cos(2*pi*tt - racro()))}
sim_square_wave       = function(tt,noise=T){rAmp()*ifelse((tt-racro())%%1<.5,1,0)+noise*rnorm(length(tt))}
sim_sqr_burst         = function(tt,noise=T){rAmp()*ifelse((tt-racro())%%1<.25,1,0)+noise*rnorm(length(tt))}

methods = list(
  sim_cosinor           = sim_cosinor,
  sim_cosinor_osc_amp   = sim_cosinor_osc_amp,
  sim_cosinor_osc_noise = sim_cosinor_osc_noise,
  sim_square_wave       = sim_square_wave,
  sim_sqr_burst         = sim_sqr_burst
)

# --------
# Panel A
# --------
tt_hd = seq(0,2,length.out=200)
signal = sim_sqr_burst(tt_hd)
data.frame(time=tt_hd,signal=signal) |> ggplot(aes(x=time,y=signal))+geom_line()


pars = expand.grid(method_name = names(methods))
set.seed(1)
dfA = lapply(seq_len(nrow(pars)), function(ii) {
  mname      = pars$method_name[ii]
  sim_method = methods[[mname]]
  Xsig       = sim_method(tt_hd,F)
  res        = data.frame(time=tt_hd,signal=Xsig,method=mname,idx=ii)
  return(res)
}
)|> rbindlist() |> data.frame()

pA=dfA |> ggplot(aes(x=time,y=signal))+geom_line()+facet_wrap(~method,scales='free')
# --------
# Panel B
# --------
pars = expand.grid(
  method_name = names(methods),
  sampling    = c( rep("random",1e2),"equispaced"),
  stringsAsFactors = FALSE
)
df = lapply(seq_len(nrow(pars)), function(ii) {
  mname      = pars$method_name[ii]
  sampling   = pars$sampling[ii]
  sim_method = methods[[mname]]
  
  if (sampling=='equispaced'){
    tt   = (0:n)/n
    tt   = tt[1:n]
  }else{
    tt = runif(n)
  }
  if (mname=='sim_cosinor_osc_noise'){ # in this case, noise should also oscillate for null
    Xnoise     = replicate(Nsamp/2,{sim_osc_noise(tt)}) |> t()
  }else{
    Xnoise     = replicate(Nsamp/2,{sim_noise(tt)}) |> t()
  }
  
  Xsig       = replicate(Nsamp/2,{sim_method(tt)}) |> t()
  X          = rbind(Xsig,Xnoise)
  resp       = c(rep(1,Nsamp/2),rep(0,Nsamp/2))
  pval       = matrixTests::row_cosinor(X,tt,period=1)$pvalue
  x          = data.frame(pval=pval,
                          resp=resp)
  rr         = pROC::roc(x,resp,pval)
  res        = data.frame(fpr=1-rr$specificities,tpr=rr$sensitivities,method=mname,sampling=sampling,idx=ii)
  return(res)
}) |> rbindlist() |> data.frame()
df |> head()
df$sampling <- factor(df$sampling, levels = c("random", "equispaced"))

method_levels <- c(
  "sim_cosinor",
  "sim_cosinor_osc_amp",
  "sim_cosinor_osc_noise",
  "sim_square_wave",
  "sim_sqr_burst"
)

dfA$method <- factor(dfA$method, levels = method_levels)
df$method  <- factor(df$method,  levels = method_levels)


pA = dfA |> ggplot(aes(x=time,y=signal))+geom_line()+facet_wrap(~method,scales='free',nrow=1)
pB = df |> ggplot(aes(x=fpr,y=tpr,group=idx,color=sampling))+geom_line()+facet_wrap(~method,nrow=1)

pA = pA + facet_wrap(~method, nrow=1,
                     labeller = labeller(method = c(
                       sim_cosinor           = "true cosinor",
                       sim_cosinor_osc_amp   = "rhythmic amplitude",
                       sim_cosinor_osc_noise = 'rhythmic noise',
                       sim_square_wave       = 'square wave',
                       sim_sqr_burst         = 'burst wave'
                     )),scales='free') +
  scale_y_continuous(labels= NULL) +
  labs(x = "time", y = "signal") +labs(color=NULL)+
  clean_theme()
pB = pB + facet_wrap(~method, nrow=1,
                     labeller = labeller(method = c(
                       sim_cosinor           = "true cosinor",
                       sim_cosinor_osc_amp   = "rhythmic amplitude",
                       sim_cosinor_osc_noise = 'rhythmic noise',
                       sim_square_wave       = 'square wave',
                       sim_sqr_burst         = 'burst wave'
                     ))) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0,.25,.5,.75,1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(x = "FPR", y = "TPR") +labs(color=NULL)+
  clean_theme()
Fig = (pA/pB)+plot_annotation(tag_levels='A')&theme(legend.position = 'bottom')
show_temp_plt(Fig,6,4)
ggsave(
  filename = "vector_figures/SuppFig04.pdf",
  plot = Fig,
  device = "pdf",
  width = 6,
  height = 3.5,
  units = "in"    
)
