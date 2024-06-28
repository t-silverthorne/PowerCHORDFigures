rowCosinor <- function(theData, zts, per=24) {
  Y <- as.matrix(theData)
  
  x1 <- sin(2*pi*zts/per)
  x2 <- cos(2*pi*zts/per)
  x0 <- rep(1, dim(Y)[2])
  X  <- cbind(x0,x1,x2)
  
  betas <- qr.solve(t(X) %*% X,t(X) %*% t(Y),tol=1e-12)
  
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