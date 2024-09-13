require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(patchwork)

am=readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
df = am@''

#df <- df %>%
#  mutate(method = if_else(method == "diffEV", method, paste0(method, lpred)))

# ncp comparisons 
df$method %>% unique()
dfgrp1 = df[,c('Nmeas','fmin','fmax','method','ncp')] %>%
  pivot_wider(names_from = method, values_from = ncp, names_prefix = "ncp_") 
p1=dfgrp1 %>% 
  ggplot(aes(y=ncp_YALMIP_bd,x=ncp_diffEV,color=fmax,shape=as.factor(Nmeas)))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_color_viridis_c()
p1


dfgrp2 = df[,c('Nmeas','fmin','fmax','method','upper')] %>%
  pivot_wider(names_from = method, values_from = upper, names_prefix = "upper_") 
dfgrp=merge(dfgrp1,dfgrp2)
p2=dfgrp %>% ggplot(aes(x=upper_YALMIP_bd,y=upper_YALMIP,color=fmax,shape=as.factor(Nmeas)))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_x_continuous(trans='log2')+ 
  scale_y_continuous(trans='log2') +scale_color_viridis_c()

p1/p2 + plot_layout(guides='collect')