# read in matlab solutions
require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
df        = read.csv2('solutions/powerCHORD_sols.csv',sep=',')
df_time   = df[,grepl('tvec',names(df))] 
df_time   = df_time%>% mutate(across(everything(),as.numeric))
df$Nmeas  = as.numeric(df$Nmeas)
df$fmin   = as.numeric(df$fmin)
df$fmax   = as.numeric(df$fmax)
df$upper  = as.numeric(df$upper)
df$lpred  = as.numeric(df$lpred)
df$MIPgap = as.numeric(df$MIPgap)
df$ncp    = as.numeric(df$ncp)
head(df)

am     = annmatrix(df_time,rann=data.frame(df[,1:7]))

df = am@''

df <- df %>%
  mutate(method = if_else(method == "diffEV", method, paste0(method, lpred)))

# ncp comparisons 
df$method %>% unique()
dfgrp1 = df[,c('Nmeas','fmin','fmax','method','ncp')] %>%
  pivot_wider(names_from = method, values_from = ncp, names_prefix = "ncp_") 

dfgrp1 %>% 
  ggplot(aes(x=ncp_YALMIP0,y=ncp_YALMIP1))+geom_point()+geom_abline(slope=1,intercept=0) 

dfgrp1 %>% 
  ggplot(aes(x=ncp_YALMIP0,y=ncp_diffEV,color=fmax))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_color_viridis_c()


# ncp comparisons 
df$method %>% unique()

dfgrp2 = df[,c('Nmeas','fmin','fmax','method','upper')] %>%
  pivot_wider(names_from = method, values_from = upper, names_prefix = "upper_") 

dfgrp2 %>% 
  ggplot(aes(x=upper_YALMIP0,y=upper_YALMIP1))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_x_continuous(trans='log2')+ 
  scale_y_continuous(trans='log2') 

dfgrp=merge(dfgrp1,dfgrp2)
dfgrp %>% ggplot(aes(x=ncp_diffEV,y=upper_YALMIP0,color=fmax))+geom_point()+geom_abline(slope=1,intercept=0)+
  scale_x_continuous(trans='log2')+ 
  scale_y_continuous(trans='log2') +facet_wrap(~Nmeas)+scale_color_viridis_c()