# setup
source('PLOSfigures/clean_theme.R')
rm()
gc()
n      = 48
tau    = c(1:n)/n -1/n
Xraw   = read.csv2('PLOSfigures/data/cutsdp_sols.csv',header = F,sep=',')
Amp    = sqrt(2) 
Nm     = 12 
tunif  = c(1:Nm)/Nm-1/Nm # equispaced
tirr   = tau[Xraw[3,]>1e-12] # good irregular

Nm1    = 8   
Nm2    = Nm-Nm1
lamb   = 4/24
t1     = c(1:Nm1)/Nm1 -1/Nm1
t2     = c(1:Nm2)/Nm2 -1/Nm2
tirr_b = c(t1*lamb,lamb +t2*(1-lamb))


Nmc     = 1e5
mt      = tunif
acrovec = 2*pi*runif(Nmc)
p24     = .5
pnull   = 0 
p4      = .5
state   = sample(x=c('circ','null','ultra'),
                 prob = c(p24,pnull,p4),
                 size=Nmc,replace=T)


# run cosinor for a list of frequencies

# plot acrophase distributions after filtering for significance as a raster plot
