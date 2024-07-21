fsize=9
theme_set(theme_classic()) 
require(dplyr)
require(matrixTests)
require(ggplot2)
require(ggplotify)
require(patchwork)

motivation_fig=function(freq,equispaced_FN,Amp=2,Nm=24,Nmc=1e4){
  tunif   = c(1:Nm)/Nm-1/Nm
  Nm1     = 2*Nm/3
  Nm2     = Nm-Nm1
  mt1     = c(1:Nm1)/Nm1-1/Nm1
  mt2     = c(1:Nm2)/Nm2-1/Nm2
  tirr    = c(mt1*.5, .5*mt2 +.5)
  acrovec = 2*pi*runif(Nmc)
  
  Yunif   = Amp*cos(outer(acrovec,2*pi*freq*tunif,'-'))+matrix(rnorm(length(tunif)*Nmc),nrow=Nmc)
  Yirr    = Amp*cos(outer(acrovec,2*pi*freq*tirr,'-'))+matrix(rnorm(length(tirr)*Nmc),nrow=Nmc)
  
  df_unif = data.frame(meas='equispaced',acro=acrovec,pval= matrixTests::row_cosinor(Yunif,tunif,1/freq) %>% {.$pvalue}) 
  df_irr  = data.frame(meas='irregular',acro=acrovec,pval= matrixTests::row_cosinor(Yirr,tirr,1/freq) %>% {.$pvalue}) 
  
  if (equispaced_FN){
    cmap_cust = c('true'=rgb(.05,0.5,.06),
                  'irregular'=rgb(.05,0.5,.06),
                  'equispaced'=rgb(.81,.06,.13))
    fn_ind  = sample(which(df_unif$pval>.05),1)
  }else{
    cmap_cust = c('true'=rgb(.05,0.5,.06),
                  'irregular'=rgb(.05,0.5,.06),
                  'equispaced'=rgb(.05,0.5,.06))
    fn_ind = sample(which(df_unif$pval<.05 & df_irr$pval<.05),1)
  }
  
  df_unif = df_unif[df_unif$pval<.05,]
  df_irr  = df_irr[df_irr$pval<.05,]
  
  df_all  = data.frame(meas='true',acro=acrovec,pval=NA)
  
  df=rbind(df_unif,df_irr,df_all)
  
  df$meas = factor(df$meas,levels=c('true','irregular','equispaced'))
  
  rad_brk = c(0,pi,2*pi)
  rad_lab = c(expression(0),
              expression(pi),
              expression(2*pi))
  plt= df %>% ggplot(aes(x=acro,fill=meas))+geom_histogram(aes(y=after_stat(density)))+facet_wrap(~meas)+
    scale_fill_manual(values=cmap_cust)
  
  plt = plt + scale_x_continuous(limits=c(0,2*pi),
                                 breaks =rad_brk,
                                 labels = rad_lab)
  plt = plt + theme(legend.position='none')
  plt = plt + labs(x=element_text('acrophase (rad)'),
                   y=element_text('density'))
  plt = plt + theme(
    strip.background=element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.25)
  )
  plt=plt+theme(text=element_text(size=fsize))
  pbot = plt
  
  
  tfine = seq(0,1,.001)
  
  acro = acrovec[fn_ind]
  dfine = data.frame(meas='true',time=tfine,signal=Amp*cos(2*pi*freq*tfine-acro))
  dunif = data.frame(meas='equispaced',time=tunif,signal=Yunif[fn_ind,])
  dirr  = data.frame(meas='irregular',time=tirr,signal=Yirr[fn_ind,])
  
  tdf = rbind(dfine,dunif,dirr)
  
  tdf$meas = factor(tdf$meas,levels=c('true','irregular','equispaced'))
  plt = tdf %>% ggplot(aes(x=time,y=signal,color=meas))+
    geom_line(data=tdf %>% filter(meas=='true'))+
    geom_point(data=tdf %>% filter(meas!='true'))+facet_wrap(~meas)+
    scale_color_manual(values=cmap_cust)+ theme(legend.position='none')
  plt = plt + labs(x=element_text('time'),
                   y=element_text('simulated expression'))
  plt = plt + scale_x_continuous(breaks=c(0,1))
  plt = plt + theme(
    strip.background=element_blank(),
    plot.margin = margin(0,0,0,0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.25)
  )
  plt=plt+theme(text=element_text(size=fsize))
  ptop = plt
  
  Fig = ptop/pbot
  return(Fig) 
}
Fig = as.ggplot(motivation_fig(1,F))/as.ggplot(motivation_fig(12,T)) + plot_annotation(tag_levels='A')

show_temp_plt=function(plt,plt_width,plt_height){
  plt_path <- tempfile(fileext = ".png")
  ggsave(plt_path, plt, width =plt_width, height = plt_height, units = "in",
         dpi = 96)
  
  viewer <- getOption("viewer")
  viewer(plt_path)
}

show_temp_plt(Fig,6,6)