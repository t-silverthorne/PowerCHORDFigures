require(tidyverse)
require(tidyr)
source('clean_figs/clean_theme.R')
Cmat = read_csv('clean_figs/data/neck_counts.csv',col_names = F)
Cmat = Cmat %>% data.frame()
Cmat = 10^Cmat
Cmat = Cmat %>% t()

rownames(Cmat)=c(4:20)
Cmat = data.frame(Cmat)
colnames(Cmat)=c('60 minutes','30 minutes','15 minutes')
Cmat$Nmeas = rownames(Cmat)

Cmat = Cmat %>% pivot_longer(cols=,names_to='grid',values_to='count')

head(Cmat)
#Cmat %>% ggplot(aes(x=))


