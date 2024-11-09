# read in diffEvolve matlab solutions
require(annmatrix)
require(dplyr)
require(ggplot2)
require(tidyr)
df        = read.csv2('diffEvolveOutput.csv',sep=',')
df_time   = df[,grepl('tvec',names(df))] 
df_time   = df_time%>% mutate(across(everything(),as.numeric))
df$Nmeas  = as.numeric(df$Nmeas)
df$fmin   = as.numeric(df$fmin)
df$fmax   = as.numeric(df$fmax)
df$upper  = as.numeric(df$upper)
df$lpred  = as.numeric(df$lpred)
df$MIPgap = as.numeric(df$MIPgap)
df$ncp    = as.numeric(df$ncp)
head(df)

am     = annmatrix(df_time,rann=data.frame(df[,1:7]))
saveRDS(am,'diffEvolveOutput.RDS')
