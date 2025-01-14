require(matrixTests)
require(parallel)
require(ggplot2)
require(data.table)
require(tidyr)
require(patchwork)
require(dplyr)
require(annmatrix)
require(devtools)
require(pROC)
require(lomb)
devtools::load_all('PowerCHORD')
pub_qual=T
mar=2

source('figures/processDiffEvolveOutput/csv_to_RDS.R')
clean_theme=function(){
  list(
    theme_classic(),
#    theme(strip.background = element_blank()),
    theme(strip.text = element_text(
      size = 7, color = "black",
      margin = margin(b = mar,r =mar,t=mar,l=mar)
    )),
    theme(
    plot.margin = margin(1,1,1,1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(vjust = 0.25)),
    theme(axis.text = element_text(size = 7)),
    theme(text=element_text(size=7))
  )
}
show_temp_plt=function(plt,plt_width,plt_height){
  if (interactive()){
    plt_path <- tempfile(fileext = ".png")
    ggsave(plt_path, plt, width =plt_width, height = plt_height, units = "in",
           dpi = 96)
    viewer <- getOption("viewer")
    viewer(plt_path)
  } 
}
rad_brk = c(0,pi/2,pi,3*pi/2,2*pi)
rad_lab = c(expression(0),
            expression(pi/2),
            expression(pi),
            expression(3*pi/2),
            expression(2*pi))

mc_cores=20 # number of cores on your machine


rowCosinor <- function(theData, zts, per=24,method=c('qr','ginv')) {
  method<-match.arg(method)
  Y <- as.matrix(theData)

  x1 <- sin(2*pi*zts/per)
  x2 <- cos(2*pi*zts/per)
  x0 <- rep(1, dim(Y)[2])
  X  <- cbind(x0,x1,x2)

  if (method=='qr'){
    betas <- qr.solve(t(X) %*% X,t(X) %*% t(Y),tol=1e-12)
  } else if (method=='ginv'){
    betas <- MASS::ginv(X)%*%t(Y) 
  }

  phases     <- atan2(betas[2,], betas[3,]) %% (2*pi)
  amplitudes <- sqrt(betas[2,]*betas[2,] + betas[3,]*betas[3,])

  fits <- t(X %*% betas)

  SStot <- rowSums((Y - rowMeans(Y))^2)
  SSres <- rowSums((fits-Y)^2)
  Rsqs  <- 1 - (SSres/SStot)

  SSmod <- SStot - SSres
  DFres <- ncol(theData) - 3
  DFmod <- 2
  MSres <- SSres / DFres
  MSmod <- SSmod / DFmod
  Fstatistic <- MSmod / MSres

  pval <- pf(Fstatistic, DFmod, DFres, lower.tail=FALSE)

  data.frame(phase=phases, amplitude=amplitudes, mesor=betas[1,],
             rsq=Rsqs, statistic=Fstatistic, pvalue=pval
  )
}
