# load in WCP designs
require(matrixTests)
require(ggplot2)
require(data.table)
require(tidyr)
N       = 18 
Amin    = 0 
Amax    = 4 
Npar    = 5e3
freq    = 3
fvec    = seq(2,4,.5) 
mt      = c(1:N)/N - 1/N
mt      =runif(N)
Amp_vec = runif(Npar,Amin,Amax)
phi_vec = runif(Npar,0,2*pi) 

eps_mat = matrix(rnorm(N*Npar),Npar,N)
Xdat    = Amp_vec*cos(outer(phi_vec,2*pi*freq*mt,'-')) + eps_mat

df_large = c(1:length(fvec)) %>% lapply(function(ii){
  freq = fvec[ii] 
  Stmat   = matrixTests::row_cosinor(Xdat,mt,period=1/freq)
  cbind(data.frame(index=c(1:Npar)),
        Stmat,
        data.frame(freq=freq,acro_true=phi_vec,Amp_true=Amp_vec))
}) %>% rbindlist() %>% data.frame()

df_plt = df_large %>% group_by(index) %>% 
  summarise(min_pval=min(pvalue),freq=freq[which.min(pvalue)],
            amp=amplitude[which.min(pvalue)],
            amp_true=ifelse(max(Amp_true)==min(Amp_true),Amp_true,NaN))

df_plt %>% filter(min_pval<.05) %>% ggplot(aes(x=amp_true,y=amp)) +
  geom_point(alpha=.45)+
  geom_abline(slope=1,intercept=0)




#df %>% ggplot(aes(x=acro_true,y=acrophase)) + geom_point()+
#  geom_abline(slope=1,intercept=0)