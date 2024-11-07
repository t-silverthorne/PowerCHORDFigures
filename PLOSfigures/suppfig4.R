source('PLOSfigures/clean_theme.R')

# load in CUTSDP solutions
n     = 48
tau   = c(1:n)/n -1/n
freqs = c(2,4,6,8,10,12)
Xraw  = read.csv2('clean_figs/data/cutsdp_sols.csv',header = F,sep=',')

###########################
# Plot of raw solutions
###########################
df=c(1:length(freqs)) %>% lapply(function(ii){
  mt = tau[as.numeric(Xraw[ii,])>0]
  mt = mt-min(mt) 
  data.frame(time=mt,fmax=freqs[ii])
}) %>% rbindlist() %>% data.frame()

dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
df$fmax_lab = paste0('{1,',df$fmax,'}')
df$fmax_lab = factor(df$fmax_lab,levels=c('{1,2}','{1,4}','{1,6}',
                                          '{1,8}','{1,10}','{1,12}'))
head(df)
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point()+
  geom_rect(data=dat_bands,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
            inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~fmax,
             ncol=1,strip.position='left')
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,2))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())

p2=plt

Fig2 = p2 

# underlying fine partition
n     = 48
tau   = c(1:n)/n -1/n
freqs = c(2,4,6,8,10,12)

# define function for testing on multifreq
check_bias=function(mt,Npar,Amin,Amax,freq_true,fvec,type,fmax_design){ # freq is true freq, fvec is test freqs
  N       = length(mt)
  Amp_vec = runif(Npar,Amin,Amax)
  phi_vec = runif(Npar,0,2*pi) 
  
  eps_mat = matrix(rnorm(N*Npar),Npar,N)
  Xdat    = Amp_vec*cos(outer(phi_vec,2*pi*freq_true*mt,'-')) + eps_mat
  
  df_large = c(1:length(fvec)) %>% lapply(function(ii){
    freq = fvec[ii] 
    Stmat   = matrixTests::row_cosinor(Xdat,mt,period=1/freq)
    Stmat$acrophase = Stmat$acrophase*2*pi*freq
    cbind(data.frame(index=c(1:Npar)),
          Stmat,
          data.frame(freq=freq,acro_true=phi_vec,Amp_true=Amp_vec))
  }) %>% rbindlist() %>% data.frame()
  
  df_out= df_large %>% group_by(index) %>% 
    summarise(min_pval=min(pvalue),freq=freq[which.min(pvalue)],
              amp=amplitude[which.min(pvalue)],
              acro=acrophase[which.min(pvalue)],
              amp_true=ifelse(max(Amp_true)==min(Amp_true),Amp_true,NaN),
              acro_true=ifelse(max(acro_true)==min(acro_true),acro_true,NaN)
    )
  
  return(cbind(df_out,data.frame(freq_true=freq_true,
                                 type=type,
                                 fmax_design=fmax_design)))
}
Npar=1e3
N=min(rowSums(Xraw)) # should be 12
mt_unif=c(1:N)/N - 1/N
mt_bad = runif(N,0,1/12) 
Amin=0
Amax=2
df=c(1:length(freqs)) %>% lapply(function(ii){
  mt = tau[as.numeric(Xraw[ii,])>1e-12]
  fvec = c(1,freqs[ii])
  rbind(check_bias(mt,Npar,Amin,Amax,1,
                   fvec,type='optimal',fmax_design=freqs[ii]),
        check_bias(mt_unif,Npar,Amin,Amax,1,
                   fvec,type='equispaced',fmax_design=freqs[ii]),
        check_bias(mt_bad,Npar,Amin,Amax,1,
                   fvec,type='fast random',fmax_design=freqs[ii]),
        check_bias(mt,Npar,Amin,Amax,freqs[ii],
                   fvec,type='optimal',fmax_design=freqs[ii]),
        check_bias(mt_bad,Npar,Amin,Amax,freqs[ii],
                   fvec,type='fast random',fmax_design=freqs[ii]),
        check_bias(mt_unif,Npar,Amin,Amax,freqs[ii],
                   fvec,type='equispaced',fmax_design=freqs[ii]))
}) %>% rbindlist() %>% data.frame()

require(ggh4x)

df$fmax_lab = paste0('{1,',df$fmax_design,'}')
df$fmax_lab = factor(df$fmax_lab,levels=c('{1,2}','{1,4}','{1,6}',
                                          '{1,8}','{1,10}','{1,12}'))

head(df)
df$signal = ifelse(df$freq_true>1,'fmax','fmin')
df$signal = factor(df$signal,levels=c('fmin','fmax'))
plt = df %>% filter(min_pval<.05)%>% 
  ggplot(aes(x=amp_true,y=amp))+
  geom_abline(slope=1,intercept = 0)+
  geom_point(alpha=.05)+
  facet_nested(fmax_lab~ type+signal)+
  ylim(c(Amin,Amax))+
  scale_y_continuous(breaks=c(Amin,Amax),
                     labels=c(Amin,Amax),
                     limits=c(Amin,Amax))+
  scale_x_continuous(breaks=c(Amin,Amax),
                     labels=c(Amin,Amax),
                     limits=c(Amin,Amax))+
  labs(x='amplitude',y='amplitude estimate')
plt=plt+clean_theme()
p1=plt

plt=df %>% filter(min_pval<.05) %>%
  ggplot(aes(x=acro_true,y=acro))+
  geom_point(alpha=.05)+
  scale_y_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  scale_x_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],labels = rad_lab[c(1,3,5)])+
  facet_nested(fmax_lab~ type+signal)+
  labs(x='acrophase (rad)',y='acrophase estimate (rad)')
plt = plt + theme(
  strip.background=element_blank(),
  plot.margin = margin(0,0,0,0),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(vjust = 0.25)
)
plt=plt+clean_theme()
p2=plt

Fig_tot = Fig2/(p1+p2)+plot_annotation(tag_levels='A')+plot_layout(heights=c(0.5,1))

show_temp_plt(Fig_tot,6,5)
ggsave(paste0('~/research/ms_powerCHORD/figures/',
              'f_cutsdp2.png'),
       Fig_tot,
       width=6,height=5,
       device='png',
       dpi=600)