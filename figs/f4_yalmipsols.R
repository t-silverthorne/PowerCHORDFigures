require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(patchwork)

am=readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
df = am@''

df <- df %>%
  mutate(method = if_else(method == "diffEV", method, paste0(method, lpred)))

# ncp comparisons 
df$method %>% unique()
dfgrp1 = df[,c('Nmeas','fmin','fmax','method','ncp')] %>%
  pivot_wider(names_from = method, values_from = ncp, names_prefix = "ncp_") 
plt=dfgrp1 %>% 
  ggplot(aes(y=ncp_YALMIP0,x=ncp_diffEV,color=fmax,shape=as.factor(Nmeas)))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_color_viridis_c(option='plasma')
plt = plt+scale_x_continuous(breaks=seq(4,24,4))+
  scale_y_continuous(breaks=seq(4,24,4))
plt = plt + labs(x='differential evolution score',
                 y='YALMIP score',
                 color='maximum frequency',
                 shape='budget')
plt=plt+guides(color=guide_colorbar(title.position='right'))
plt=plt+guides(shape=guide_legend(title.position='right'))
plt=plt+theme(text=element_text(size=fsize),legend.direction='vertical',
                legend.title = element_text(angle = 90,hjust=0.5))
p1=plt




dfgrp2 = df[,c('Nmeas','fmin','fmax','method','upper')] %>%
  pivot_wider(names_from = method, values_from = upper, names_prefix = "upper_") 
dfgrp=merge(dfgrp1,dfgrp2)
plt=dfgrp %>% ggplot(aes(x=ncp_diffEV,y=upper_YALMIP0,color=fmax,shape=as.factor(Nmeas)))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_x_continuous(trans='log2')+ 
  scale_y_continuous(trans='log2') +scale_color_viridis_c(option='plasma')
plt = plt+scale_y_continuous(breaks=seq(8,50,8))+
  scale_x_continuous(breaks=seq(4,24,4))
plt = plt + labs(x='differential evolution score',
                 y='YALMIP upper bound',
                 color='maximum frequency',
                 shape='budget')
plt=plt+guides(color=guide_colorbar(title.position='right'))
plt=plt+guides(shape=guide_legend(title.position='right'))
plt=plt+theme(text=element_text(size=fsize),legend.direction='vertical',
                legend.title = element_text(angle = 90,hjust=0.5))
plt

p2=plt

Fig = p1/p2 + plot_layout(guides='collect') + plot_annotation(tag_levels='A')
show_temp_plt(Fig,6,3.5)


ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f4_yalmipsols.png'),
       Fig,
       width=6,height=3.5,
       device='png',
       dpi=600)