require(matrixTests)
require(gcrma)
require(GeneCycle)



# load data
df = read.csv('data_analysis/GSE11923_series_matrix.txt',
              skip='65',header=F,sep='\t',row.names=1)
names(df)=paste0('CT',18:65)
df=df[!is.na(rowSds(as.matrix(df))>0),]
head(df)

# naive cosinor
cs24   = matrixTests::row_cosinor(df,c(18:65),period = 24) 
q_cs24 = p.adjust(cs24$pvalue,method='fdr')
cs12   = matrixTests::row_cosinor(df,c(18:65),period = 12) 
q_cs12 = p.adjust(cs12$pvalue,method='fdr')
cs8    = matrixTests::row_cosinor(df,c(18:65),period = 8) 
q_cs8  = p.adjust(cs8$pvalue,method='fdr')

# Fisher g test 
fg      = GeneCycle::fisher.g.test(t(as.matrix(df)))
q_gtest = p.adjust(fg,method='fdr')

# find genes of each periodicity
((q_gtest<.05) & (q_cs24<.05)) |> summary()
((q_gtest<.05) & (q_cs12<.05)) |> summary()
((q_gtest<.05) & (q_cs8<.05))  |> summary()

# check overlap with those reported in original paper
orig_24 = read.csv2('data_analysis/pgen.1000442.s012.csv',sep=',')
orig_12 = read.csv2('data_analysis/pgen.1000442.s013.csv',sep=',')
orig_8  = read.csv2('data_analysis/pgen.1000442.s014.csv',sep=',')


( rownames(df)[((q_gtest<.05) & (q_cs24<.05))] %in% orig_24[,1] ) |> summary() 
( rownames(df)[((q_gtest<.05) & (q_cs12<.05))] %in% orig_12[,1] ) |> summary()
( rownames(df)[((q_gtest<.05) & (q_cs8<.05))]  %in% orig_8[,1]  ) |> summary()

( orig_24[,1]  %in% rownames(df)[((q_gtest<.05) & (q_cs24<.05))] ) |> summary() 
( orig_12[,1]  %in% rownames(df)[((q_gtest<.05) & (q_cs12<.05))] ) |> summary()
( orig_8[,1]   %in% rownames(df)[((q_gtest<.05) & (q_cs8<.05))]  ) |> summary()

fisher.test()





# COSOPT
# pull this from the cosopt repo and install if needed
# install.packages("./data_analysis/cosopt_0.3.1.tar.gz",repos = NULL, type = "source")
#co = cosopt(data=as.numeric(df[1,]),
#            sigma=rep(0.1,48),
#            timepoints = c(18:65))
