evalMinEig <- function(t,freq){
  A      = matrix(c(0,0,1,0,0,1),nrow=3,byrow=T)
  Xr     = matrix(c(cos(2*pi*freq*t),sin(2*pi*freq*t)),ncol=2)
  D      = t(Xr)%*%Xr
  b      = matrix(c(sum(cos(2*pi*freq*t)),sum(sin(2*pi*freq*t))),ncol=1)
  invB   = D - b%*%t(b)/length(t)
  return(eigen(invB) |> (\(x) x$values)() |> min())
}