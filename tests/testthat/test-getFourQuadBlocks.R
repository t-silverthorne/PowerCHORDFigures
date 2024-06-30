test_that("symmetry of Cm12 Cm21", {
  freq  = runif(1)*10
  Nfine = sample(c(100:300),1)
  Nmeas = sample(c(3:50),1)
  mats  = getFourQuadBlocks(freq,Nfine,Nmeas)
  expect_equal(max(abs(mats$Cm12-mats$Cm21)),0)
})

#TODO: check agrees with B matrix for different freqs
