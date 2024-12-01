# Import dataset
require(matrixTests)
require(matrixStats)
require(ggplot2)
require(dplyr)
require(data.table)
require(GeneCycle)
require(devtools)
require(ggh4x)
require(dplyr)
require(patchwork)
devtools::load_all("PowerCHORD")
# load data 
df = read.csv('data_analysis/GSE11923_series_matrix.txt',
              skip='65',header=F,sep='\t',row.names=1)
names(df)=paste0('CT',18:65)
df=df[!is.na(rowSds(as.matrix(df))>0),]
head(df)

####################
# Get FDR statistic on sub-sampling
####################

# measurement schedules
mt_full  = c(18:65)
dt       = 4
mt_sub   = seq(18,65,dt)

fname=paste0('data_analysis/mu_opt',48/dt,'.txt')
mu_opt   = as.numeric(t(read.csv(fname,header=F))) 
mt_opt   = mt_full[mu_opt>0]

length(mt_sub)
length(mt_opt)


pars = expand.grid(sampling = c('full','sub','optimal'),
                   per = c(8,12,24))

sdf = c(1:dim(pars)[1]) |> lapply(function(ind){
  samp = pars[ind,]$sampling 
  per  = pars[ind,]$per
  
  if (samp == 'full'){
    dfloc=df
    mtloc=c(18:65)
  }else if (samp == 'sub'){
    dfloc=df[,mt_full %in% mt_sub] 
    mtloc=mt_sub
  }else if (samp == 'optimal'){
    dfloc=df[,mt_full %in% mt_opt] 
    mtloc=mt_opt
  }
  cs   =  matrixTests::row_cosinor(dfloc,mtloc,period = per) 
  qval = p.adjust(cs$pvalue,method='fdr')
  pdf = data.frame(probe   = rownames(dfloc),
                   p_osc   = cs$pvalue,
                   q_osc   = qval,
                   acro = cs$acrophase,
                   amp = cs$amplitude,
                   per = per,
                   sampling=samp
  )
}
) |> rbindlist() |> data.frame()

sdf_copy = sdf
source('PLOSfigures/clean_theme.R')
sdf$samp_wrap = factor(sdf$sampling, levels = c("full", "optimal", "sub"),
                       labels = c("full (N=48)", "optimal (N=12)", 
                                  "sub-sampled (N=12)"))
sdf$per_wrap = factor(sdf$per, levels = c(8,12,24),
                      labels = c("T=8 hr", "T=12 hr", "T=24 hr"))

sdf$acro[(sdf$per==8) & (sdf$sampling=='sub')]=NA
p1 = sdf |> filter(p_osc<.05) |>  ggplot(aes(x=2*pi*acro/per))+
  geom_histogram()+facet_grid(per_wrap~samp_wrap,scales = 'free')+
  scale_y_continuous(trans='log2')+clean_theme()+
  labs(x='acrophase (rad)',y='log2(count)')

tdf = rbind(data.frame(time = mt_full,sampling='full'),
        data.frame(time = mt_sub,sampling='sub'),
        data.frame(time = mt_opt,sampling='optimal'))

tdf$sampfac = factor(tdf$sampling, levels = c("full", "optimal", "sub"),
                       labels = c("full (N=48)", "optimal (N=12)", 
                                  "sub-sampled (N=12)"))
plt = tdf |> ggplot(aes(x=time,y=1))+
  geom_point(size=.1)+facet_wrap(~sampfac)+labs(x='time (hr)',y=NULL)+
  clean_theme()
plt = plt + clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
p2=plt

Fig = (p2/p1) + plot_layout(heights = c(1,3))+plot_annotation(tag_levels='A')
Fig
ggsave(paste0('PLOSfigures/',
              'fig7.png'),
       Fig,
       width=6,height=4,
       device='png',
       dpi=600)