Nsamp=3
Xmat = data.frame(pvalue=runif(Nsamp*2),
                  freq=c(rep(1,Nsamp),rep(2,Nsamp)),
                  idx=c(c(1:Nsamp),c(1:Nsamp)))
Xmat
Xmat %>% group_by(idx) %>% 
  summarise(min_pval=min(pvalue),freq=freq[which.min(pvalue)])


  mutate(row=row_number()) %>% 
  pivot_wider(names_from=idx,values_from=c('freq','pvalue')) %>% 
Xmat %>% pivot_wider(names_from = idx,values_from=c('freq','pvalue')) %>%summarise(min()) 
Xmat_wide <- Xmat %>%
  pivot_wider(names_from = idx, values_from = c('freq', 'pvalue'))

