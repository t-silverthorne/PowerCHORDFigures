require(devtools)
require(gurobi)
require(data.table)
load_all()

Nfine=144

pars=expand.grid(MIPFocus=c(2),MIQCPMethod=c(0,1))

df=c(1:dim(pars)[1]) %>% lapply(function(ind){
res=runPcNCQ(freq=10,Nfine=Nfine,Nmeas=48,con_mode='exact',WorkLimit=100,
             MIPFocus=pars[ind,]$MIPFocus,MIQCPMethod=pars[ind,]$MIQCPMethod)
return(cbind(pars[ind,],data.frame(gap=res$sol$mipgap*100)))
}) %>% rbindlist() %>% data.frame()

df[which.min(df$gap),]

#sol = res$sol
#
#a=sol$x[Nfine+1]
#b=sol$x[Nfine+2]
#c=sol$x[Nfine+3]
#matrix(c(a,b,b,c),nrow=2) %>% eigen() %>% {.$values} %>% min()
#sol$x[Nfine+8]*sol$x[Nfine+10]

