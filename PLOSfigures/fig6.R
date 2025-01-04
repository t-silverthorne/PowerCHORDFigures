source('PLOSfigures/clean_theme.R')
require(ggrepel)

# load in CUTSDP solutions
N     = 12
n     = 48
tau   = c(1:n)/n -1/n
freqs = c(2,4,6,8,10,12)
Xraw  = read.csv2('PLOSfigures/data/cutsdp_sols.csv',header = F,sep=',')

###########################
# Plot of raw solutions
###########################
topt  = tau[as.numeric(Xraw[6,])>0]
tunif = c(1:N)/N - 1/N
length(topt)==length(tunif)

df = rbind(data.frame(time=topt,type='optimal'),
      data.frame(time=tunif,type='equispaced'))

dat_bands = data.frame(start=c(0:12)*2,end=(c(0:12)*2+1))
plt = df %>% ggplot(aes(x=24*time,y=1))+geom_point()+
  geom_rect(data=dat_bands,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf),alpha=.4,
            inherit.aes = F,fill=c('lightblue'))+
  facet_wrap(~type,
             ncol=1,strip.position='top')
plt=plt+clean_theme()
plt = plt + theme(axis.title.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + theme(axis.text.y=element_blank()) 
plt = plt + labs(x=element_text('time (hr)'))
plt = plt+scale_x_continuous(limits=c(0,24),breaks=seq(0,24,2))
plt = plt+theme(axis.line.y = element_blank())
plt = plt+theme(axis.ticks.y = element_blank())
p1=plt

###########################
# Cycle diagram
###########################
plt=data.frame(time=topt,rr=2*(1+((24*topt)%%2))) |> 
  ggplot(aes(x=time*2*pi,y=rr))+
  geom_point(size=1)+coord_polar(theta='x',clip='off')+theme_minimal()+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
                     labels = c('0','6','12','18'),
                     limits = c(0,2*pi))+
  geom_hline(yintercept = seq(2, 5, by = 1), linetype = "dashed", 
             color = "black") +
  scale_y_continuous(labels=c(),limits=c(1,6),breaks=NULL)+
  labs(x='phase (along 24hr cycle)',y=NULL)
plt=plt+clean_theme()
plt = plt + theme(
  plot.margin = margin(0,0,0,0),
  #axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line = element_blank()
)
plt
q1=plt

plt=data.frame(time=2*pi*((24*topt)%%2)/2,phase=24*topt+1) |> 
  ggplot(aes(x=time,y=phase))+
  geom_point(size=1)+coord_polar(theta='x')+theme_minimal()+
  scale_x_continuous(breaks=c(0,pi/2,pi,3*pi/2),
                     labels = c('0','1/2','1','3/2'),
                     limits = c(0,2*pi))+
  geom_hline(yintercept = seq(4, 24, by = 6), 
             linetype = "dashed", color = "black") +
  scale_y_continuous(labels=c(),limits=c(0,24),breaks=NULL)+
  labs(x='phase (along 2hr cycle)',y=NULL)
plt=plt+clean_theme()
plt = plt + theme(
  plot.margin = margin(0,0,0,0),
  #axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line = element_blank()
)
plt
q2=plt

###########################
# ncp across freqs 
###########################
fvec = seq(0,16,.01)
pars = expand.grid(freq=fvec,type=c('optimal','equispaced'))

edf = c(1:dim(pars)[1]) %>% lapply(function(ii){
  x=pars[ii,]
  freq=x[['freq']]
  if(x[['type']]=='equispaced'){
    mt = tunif
  }else if(x[['type']]=='optimal'){
    mt = topt 
  }else{
    stop('unknown type')
  }
  return(cbind(pars[ii,],data.frame(emin =evalMinEig(mt,freq))))
}) %>% rbindlist() %>% data.frame()

color_scale= c('optimal'='#F8766D',
               'equispaced'='#00BFC4')
plt = edf %>% ggplot(aes(x=freq,y=emin,group=type,color=type))+geom_line()+
  geom_vline(xintercept = 12,linetype='dashed')+
  geom_vline(xintercept = 1,linetype='dashed')+
  scale_x_continuous(limits=c(0,16),breaks=seq(0,16,4))+
  scale_color_manual(values=color_scale)+
  labs(x='frequency (cycles/day)',y='noncentrality')+
  guides(color=guide_legend(title=NULL))
plt = plt+clean_theme()
plt = plt+theme(legend.position='bottom')
p2  = plt

data.frame(time = (topt %% (1/12))) %>% ggplot(aes(x=time,y=1))+geom_point()


###########################
# Pareto frontier  
###########################

nrep = 1e4

#random designs
rdf=c(1:nrep) %>% lapply(function(ii){
  mt = runif(N)  
  return(data.frame(type='random',
                    eig1=evalMinEig(mt,1),
                    eig12=evalMinEig(mt,12)))
}) %>% rbindlist() %>% data.frame()
# equispaced
rdf=rbind(rdf,data.frame(type='equispaced',eig1=evalMinEig(tunif,1),eig12=evalMinEig(tunif,12)),
data.frame(type='optimal',eig1=evalMinEig(topt,1),eig12=evalMinEig(topt,12)))

alpha_scale = c('optimal'=1,'equispaced'=1,'random'=.1)
size_scale  = c('optimal'=2.5,'equispaced'=2.5,'random'=.7)
type_scale = c('optimal'=18,'equispaced'=18,'random'=19)
color_scale= c('optimal'='#F8766D',
               'equispaced'='#00BFC4',
               'random'='black')
plt = rdf %>% ggplot(aes(x=eig1,y=eig12,color=type,size=type,
                       alpha=type,shape=type))+
  geom_point()+
  scale_size_manual(values=size_scale)+
  scale_alpha_manual(values = alpha_scale)+
  scale_color_manual(values=color_scale)+
  scale_shape_manual(values=type_scale)+
  geom_hline(yintercept = 6,color='purple',linetype='dashed')+
  geom_vline(xintercept = 6,color='purple',linetype='dashed')+
  geom_label_repel(data = subset(rdf, type %in% c("equispaced", "optimal")),
                  aes(label = type),show.legend=F)+
  labs(x='noncentrality (f=1)',y='noncentrality (f=12)')+
  guides(color=guide_legend(title=NULL),
         size=guide_legend(title=NULL),
         alpha=guide_legend(title=NULL),
         shape=guide_legend(title=NULL)
         )
plt = plt+clean_theme()
plt = plt+theme(legend.position='bottom')
p3  = plt
p3


##############
# Calendar
##############
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggplotify)
# Generate data for a generic 4x7 calendar grid
calendar_data <- expand.grid(Week = 1:4, Day = 1:7) %>%
  mutate(DayType = ifelse(Day == 6, "Saturday", "Other"))

# Define positions for the dots on a circle
saturday_dots <- calendar_data %>%
  filter(DayType == "Saturday") %>%
  tidyr::expand_grid(
    Angle = seq(0, 2 * pi, length.out = 5)[-5]
  ) %>%
  mutate(
    radius = 0.3, 
    DotX = Day + radius * cos(Angle),
    DotY = Week + radius * sin(Angle)
  )

# Define the circle outline
circle_outline <- calendar_data %>%
  filter(DayType == "Saturday") %>%
  tidyr::expand_grid(
    Angle = seq(0, 2 * pi, length.out = 100)
  ) %>%
  mutate(
    radius = 0.3, 
    CircleX = Day + radius * cos(Angle),
    CircleY = Week + radius * sin(Angle),
    Group = interaction(Day, Week)
  )

# Function to create a calendar plot
create_calendar <- function(data_grid, dots_data = NULL, circles_data = NULL) {
  p <- ggplot(data_grid, aes(x = Day, y = Week)) +
    geom_tile(color = "black", fill = "white") +
    coord_fixed() +
    scale_y_reverse() + # Reverse y-axis
    scale_x_continuous(breaks = 1:7) + # X-axis labels for each day (1 to 7)
    labs(x = "Day", y = "Week") + # Add axis titles
    theme_minimal() +
    theme(
      axis.title = element_text(), # Axis titles
      plot.title = element_blank(), # Remove panel titles
      axis.text.x = element_text(), # Size of x-axis labels
      axis.text.y = element_text()  # Size of y-axis labels
    )
  
  # Add dots if provided
  if (!is.null(dots_data)) {
    p <- p + 
      geom_point(size=.2,data = dots_data, aes(x = DotX, y = DotY), color = "black", size = 2)
  }
  
  # Add circles if provided
  if (!is.null(circles_data)) {
    p <- p +
      geom_path(
        linewidth=.2,
        data = circles_data, 
        aes(x = CircleX, y = CircleY, group = Group),
        color = "black",
      )
  }
  
  p=p+clean_theme()
  p
}

# First row of calendars
calendar_with_dots <- create_calendar(
  calendar_data,
  dots_data = saturday_dots,
  circles_data = circle_outline
)

empty_calendar <- create_calendar(calendar_data)

first_row <- (calendar_with_dots / empty_calendar / empty_calendar)
first_col <- (calendar_with_dots | empty_calendar | empty_calendar)

# Second row of calendars, dots separated by weeks
dots_first_week <- create_calendar(
  calendar_data,
  dots_data = saturday_dots %>% filter(Week == 1),
  circles_data = circle_outline %>% filter(Week == 1)
)

dots_second_week <- create_calendar(
  calendar_data,
  dots_data = saturday_dots %>% filter(Week == 2),
  circles_data = circle_outline %>% filter(Week == 2)
)

dots_third_fourth_weeks <- create_calendar(
  calendar_data,
  dots_data = saturday_dots %>% filter(Week %in% c(3, 4)),
  circles_data = circle_outline %>% filter(Week %in% c(3, 4))
)

second_row <- (dots_first_week / dots_second_week / dots_third_fourth_weeks) 
second_col <- (dots_first_week | dots_second_week | dots_third_fourth_weeks) 
frgp = as.ggplot(first_row)
srgp = as.ggplot(second_row)

# Combine the rows
full_cal  = ( frgp / srgp)
####################################

Fig=(p1|p2)/ ( ((q1/q2) |(first_col/second_col)) +plot_layout(widths=c(1.5,3))) + 
  plot_layout(heights=c(1.5,8))+
  plot_annotation(tag_levels=list(c('A','B','C','D','E','','','F')))
show_temp_plt(Fig,6,4.5)
ggsave('PLOSfigures/fig6.png',
       Fig,
       width=6,height=4.5,
       device='png',
       dpi=600)

####################################
# Supp fig part
####################################

#logscale freq to show circadian/circalunar/circannual spectrum

pers = c(24, 24*7*4, 24*7*4*12)

# Log-transform the values
log_pers <- log10(pers)

# Create a sequence of evenly spaced points in the log space
log_interp <- seq(from = log_pers[1], to = log_pers[3], length.out = 1e3)

year_hr = pers[3]

mt0 = c(0,6,12,18)
mt = c(mt0,mt0+7*24,mt0+14*24,mt0+21*24)
mt0 = mt
mt  = c(mt0,mt0+3*4*7*24,mt0+6*4*7*24,mt0+9*4*7*24)
mt  = mt/year_hr


freq1 = 1
freq2 = year_hr/24/7/4 
freq3 = year_hr/24

Nacro     = 2^8+1
acros     = seq(0,2*pi,length.out=Nacro)
acros     = acros[1:(length(acros)-1)]


pars = expand.grid(freq_type  = c('circadian','circalunar','circannual'),
                   freq_delta = seq(0.95,1.05,length.out=1e3),
                   acro       = acros)

pwr_vec =c(1:dim(pars)[1]) |> mclapply(mc.cores=mc_cores,function(ii){
  x          = pars[ii,]
  acro       = x$acro
  freq_type  = x$freq_type
  freq_delta = x$freq_delta 
  freq=NaN
  if (freq_type=='circadian'){
    freq=year_hr/24
  }else if (freq_type=='circalunar'){
    freq=year_hr/24/7/4
  }else if (freq_type=='circannual'){
    freq=1
  }
  freq=freq*freq_delta
  return(evalPower(mt,acro=acro,freq=freq,Amp=1/sqrt(2)) )
})
pars$power = as.numeric(pwr_vec)

plt = pars |> ggplot(aes(y=acro,x=freq_delta,fill=power))+geom_raster()+
  scale_fill_viridis_c()+
  facet_wrap(~freq_type,ncol=1)+
  scale_y_continuous(limits=c(0,2*pi),breaks =rad_brk[c(1,3,5)],
                     labels = rad_lab[c(1,3,5)])+
  scale_x_continuous(limits=c(0.95,1.05),breaks=c(0.95,1,1.05))+
  labs(y='acrophase (rad)',x='relative frequency')+clean_theme()
plt = plt + theme(legend.direction='horizontal',legend.position='bottom')
r1=plt

length(mt)
evalMinEig(mt,year_hr/24)
evalMinEig(mt,year_hr/24/7/4)

Figsup = p3+r1 + plot_annotation(tag_levels='A')+
  plot_layout(widths = c(1,2))
show_temp_plt(Figsup,6,3.5)

ggsave('PLOSfigures/fig6sup.png',
       Figsup,
       width=6,height=3.5,
       device='png',
       dpi=600)

