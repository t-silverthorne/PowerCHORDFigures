# Composite figure: peerReviewFig1 + Panel A from freePeriod_multipanel
require(ggplot2)
require(dplyr)
require(patchwork)
source('figures/clean_theme.R')

# -----------------------------
# Build Fig (from peerReviewFig1.R)
# -----------------------------
df1        = read.csv('additional_figures/data/prF1c_n48_modereal.csv',F)
names(df1) = c('Amp','fmax','freq','power')
df1        = cbind(df1,data.frame(type='equispaced'))
df2        = read.csv('additional_figures/data/prF1c_cheb_n48_modereal.csv',F)
names(df2) = c('Amp','fmax','freq','power')
df2        = cbind(df2,data.frame(type='irregular'))
df3        = read.csv('additional_figures/data/prF1c_diffEv_n48_modereal.csv',F)
names(df3) = c('Amp','fmax','freq','power')
df3        = cbind(df3,data.frame(type='diffEv'))
df         = rbind(df1,df2,df3)
df$type    = factor(df$type, 
                    levels = c('equispaced','irregular','diffEv'),
                    labels=c('equispaced','permutation bound','fixed-period heuristic'))
df$fmax = factor(df$fmax, levels = c(12,16,24),
            labels = c(expression(f[max]==N/4),
                       expression(f[max]==N/3),
                      expression(f[max]==N/2)))
df  = df |> filter(Amp>1)
df$Amp  = factor(df$Amp, levels = c(1.5,2),
            labels = c(expression(Amp==1.5),
                       expression(Amp==2)))

Fig = df |> ggplot(aes(x=freq,y=power,color=type,group=type))+
  facet_grid(Amp~fmax, scales='free_x',labeller=label_parsed)+
  geom_line(linewidth=0.5)+geom_point(size=0.5)+labs(x='frequency')+
  coord_cartesian(xlim = c(1, NA))+clean_theme()+
  labs(color = NULL)+
  scale_color_manual(values = c(
    'equispaced' = '#619CFF',  # ggplot default blue
    'permutation bound' = '#F8766D',  # ggplot default red
    'fixed-period heuristic' = '#00BA38'   # ggplot default green
  ))+
  scale_x_continuous(
    limits = c(1, NA),
  breaks = function(x) {
    # x is the range of the axis, x[2] is the max
    max_x = ceiling(max(x, na.rm=TRUE))
    breaks = c(1,seq(4, max_x, by=4))
  })& theme(legend.position='bottom')

# ---------------------------------------------
# Build Panel A (from freePeriod_multipanel.R)
# ---------------------------------------------
dfa  = read.csv('additional_figures/data/results_prF1a_equi_1.csv',header=F)
dfa2 = read.csv('additional_figures/data/results_prF1a_equi_2.csv',header=F)
dfa3 = read.csv('additional_figures/data/results_prF1a_equi_3.csv',header=F)
df_a = rbind(dfa,dfa2,dfa3)
names(df_a) = c('Nmeas','fmin','fmax','Amp','acro','freq','pfree','pfixed')
df_a$fmax = factor(df_a$fmax, levels = c(3,4,6),
                 labels = c(expression(f[max]==N/4),
                            expression(f[max]==N/3),
                            expression(f[max]==N/2)))
p_a = df_a |> ggplot(aes(x=pfixed,y=pfree,color=freq))+geom_point(size=.3)+
  facet_wrap(~fmax,labeller=label_parsed)+  scale_color_viridis_c(limits=c(1,6),
                        name   = "frequency",
                        )+clean_theme()+labs(x='fixed period power',y='free period power')+
  scale_x_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1),
    labels = c("0", "0.25", "0.5", "0.75", "1"))+
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1),
    labels = c("0", "0.25", "0.5", "0.75", "1"))
p_a = p_a + guides(color='none')

# -----------------------------
# Combine and save
# -----------------------------
Combined = (Fig/p_a) + plot_annotation(tag_levels='A')+
  plot_layout(heights=c(1.5,1),guides='collect') & theme(legend.position='bottom')
ggsave(
  filename = "vector_figures/Main03.pdf",
  plot = Combined,
  device = "pdf",
  width = 6,
  height = 4,
  units = "in"    
)