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

Cmat$Nmeas = as.numeric(Cmat$Nmeas)
Cmat = Cmat %>% pivot_longer(cols=c('60 minutes','30 minutes','15 minutes'),
                             names_to='grid',values_to='count')

head(Cmat)
plt=Cmat %>% ggplot(aes(x=Nmeas,y=count,group=grid,
                    color=grid))+
  geom_line()+geom_point()+
scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
              labels = scales::trans_format("log10", scales::math_format(10^.x)))+
scale_x_continuous(breaks=seq(4,20,4))+
  labs(x='sample size',y='number of design\nequivalence classes',
       color='grid spacing')
plt=plt+clean_theme()
show_temp_plt(plt,6,2)

ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'fig_necklace.png'),
       plt,
       width=6,height=2,
       device='png',
       dpi=600)