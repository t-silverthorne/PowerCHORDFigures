require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
require(patchwork)


am=readRDS('figures/sec3p2_data/powerCHORD_even_sols.RDS')
df = am@''

dfgrp1 = df[,c('Nmeas','fmin','fmax','method','ncp')] %>%
  pivot_wider(names_from = method, values_from = ncp, names_prefix = "ncp_") 
dfgrp2 = df[,c('Nmeas','fmin','fmax','method','upper')] %>%
  pivot_wider(names_from = method, values_from = upper, names_prefix = "upper_") 
dfgrp=merge(dfgrp1,dfgrp2)

dfgrp[which(dfgrp$ncp_diffEVCR>dfgrp$upper_YALMIP_bd ),]
am@'' %>% names()
am@method %>% unique()
uu=am[am@Nmeas==16 & am@fmin ==6 & am@fmax==20 & am@method=='diffEVCR' ,]
mt = as.numeric(uu)[!is.nan(uu)]
fvec = seq(6,20,.5)

fvec %>% sapply(function(freq){
  evalMinEig(mt,freq)
}) %>% min()


summary(dfgrp$ncp_diffEVCR>dfgrp$Nmeas/2)
summary(dfgrp$Nmeas/2>dfgrp$upper_YALMIP_bd )