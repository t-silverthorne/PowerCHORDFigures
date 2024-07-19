require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(patchwork)

am=readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
df = am@''
df$method %>% unique() |> print()

df = df[df$method%in%c('diffEV','diffEVCR'),]

dfgrp1 = df[,c('Nmeas','fmin','fmax','method','ncp')] %>%
  pivot_wider(names_from = method, values_from = ncp, names_prefix = "ncp_") 
p1=dfgrp1 %>% 
  ggplot(aes(y=ncp_diffEV,x=ncp_diffEVCR,color=fmax,shape=as.factor(Nmeas)))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_color_viridis_c()