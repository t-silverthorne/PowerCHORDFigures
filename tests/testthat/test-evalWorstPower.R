test_that("phase discretization matches eigenvalue method", {
  param = list(Amp=runif(1)*2,freq=runif(1,0.5,12),acro=NaN)
  Nmeas = sample(6:12,1)
  mt    = c(1:Nmeas)/Nmeas-1/Nmeas
  mp1=evalWorstPower(mt,param,method='test')
  mp2=evalWorstPower(mt,param,method='eig')
  expect_equal(mp1,mp2,tolerance =1e-6)
})
