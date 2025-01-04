require(data.table)
require(ggplot2)
require(dplyr)
require(devtools)
devtools::load_all('PowerCHORD')
data.frame(mt=(tirr %% (1/48))) |> ggplot(aes(x=mt,y=1))+geom_point()

(24*tirr) %% 4
(24*tirr) %% 8  

( (48*tirr) %% 6 )


# known example
f1=1 # T1 = 24
f2=6 # T2 = 4
tbase = c(0,8,16)
mt = c(tbase,tbase+1,tbase+2,tbase+3)/24
evalMinEig(mt,1)
evalMinEig(mt,6)

acrovec=seq(0,2*pi,length.out=2^10)

# new example
f1 = 1 # T1 = 24
f2 = 4 # T2 = 6
tbase=6*c(0:4)/4
tbase=tbase[1:4]
mt = c(tbase,tbase+6,tbase+12,tbase+18)/24
evalMinEig(mt,1)
evalMinEig(mt,4)

pwr = acrovec |> lapply(function(acro){
  evalPower(mt,freq=4,acro,Amp=1)
}) |> unlist()
data.frame(acro=acrovec,power=pwr) |> ggplot(aes(x=acro,y=pwr))+geom_line()
pwr = acrovec |> lapply(function(acro){
  evalPower(mt,freq=1,acro,Amp=1)
}) |> unlist()
data.frame(acro=acrovec,power=pwr) |> ggplot(aes(x=acro,y=pwr))+geom_line()

# harder example
T1=7*8
T2=4
n1=7
k1=2
n2=4

tbase = T2*c(0:n2)/n2
tbase=tbase[1:n2]

mt = tbase
for (ii in 1:((k1*n1)-1)){
  mt=c(mt,tbase+ii*k1*T2)
}
length(mt)
evalMinEig(mt/T1,1)
evalMinEig(mt/T1,T1/T2)

# recursive 2 step
nf = 2
Tvec = c(4,7*8)
kvec = c(2)
nvec = c(4,7)

get_meas_times=function(kk,Tvec,nvec,kvec,mt=NULL){
  if (kk==1){
    tbase = Tvec[kk]*c(0:nvec[kk])/nvec[1]
    tbase = tbase[1:nvec[kk]]
    mt    = tbase
  }else{
    tbase=mt
  } 
  
  for (ii in 1:(nvec[kk+1]-1)){
    mt = c(mt,tbase+ii*kvec[kk]*Tvec[kk]) 
  }
  if (kk<length(Tvec)-1){
   get_meas_times(kk+1,Tvec,nvec,kvec,mt) 
  }else{
    return(mt)
  } 
}
check_meas_times=function(mt,Tvec){
  nf =length(Tvec)
  c(1:nf) |> lapply(function(kk){
    data.frame(period=Tvec[kk],min_eig = evalMinEig(mt,Tvec[nf]/Tvec[kk])) 
  }) |> rbindlist() |> data.frame()
}
mt = get_meas_times(1,Tvec,nvec,kvec,NULL)
mt = mt/Tvec[nf]
check_meas_times(mt,Tvec)





# recursive 2 step
Tvec = c(1,8,24)
nf   = length(Tvec)
nvec = c(4,4,3)
kvec = rep(NaN,nf-1)
for (ii in c(1:nf-1)){
  kvec[ii] = Tvec[ii+1]/nvec[ii+1]/Tvec[ii]
}
mt = get_meas_times(1,Tvec,nvec,kvec,NULL)
mt = mt/Tvec[nf]
data.frame(time=mt) |> ggplot(aes(x=time ,y=1))+geom_point()
check_meas_times(mt,Tvec)


#Tvec = c(24,24*7*4,24*7*4*12)
#nvec_best = c(10,5,9)
Tvec = c(24,24*7*4,24*7*4*12)
nvec = c(4,4,6)
nf   = length(Tvec)
kvec = rep(NaN,nf-1)
for (ii in c(1:nf-1)){
  kvec[ii] = Tvec[ii+1]/nvec[ii+1]/Tvec[1]
}
mt = get_meas_times(1,Tvec,nvec,kvec,NULL)
mt = mt/Tvec[nf]
mt = mt
max(mt)
data.frame(time=mt%%1) |> 
  ggplot(aes(x=(time)*12 ,y=1))+geom_point(size=.1)
check_meas_times(mt,Tvec) 
check_meas_times(mt%%1,Tvec) 
length(mt)
length(unique(mt))