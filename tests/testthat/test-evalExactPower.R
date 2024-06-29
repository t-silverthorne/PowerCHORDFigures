require(matrixTests)
test_that("compare with monte carlo", {
  param = list(freq=runif(1,0.5,3),
               Amp=1+.1*runif(1),
               acro=2*pi*runif(1))
  mt = c(1:25)/25-1/25
  Nmc   = 1e4
 
  # optimal choice of Nperm based on Boos Zhang heuristic 
  malpha  = 10
  al_val  = 1/malpha
  Npcands = malpha*c(1:1000)-1
  Nperm   = Npcands[which.min(abs(Npcands-8*sqrt(Nmc)))]
  
  pwr_exact = evalExactPower(mt,param,al_val)
  pwr_MC1 = evalMonteCarloPower(mt,param,Nmc,al_val,method='Ftest')
  expect_equal(pwr_exact,pwr_MC1,tolerance = 1e-2)
})