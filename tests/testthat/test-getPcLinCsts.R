test_that("Measurement budget constraint", {
  Nmeas=sample(c(12:48),1)
  Nfine=120

  LC=getPcLinCsts(Nmeas,Nfine,Inf)
  expect_equal(LC$A,matrix(c(rep(1,Nfine),0),nrow=1))
  expect_equal(LC$sense_list,list('='))
  expect_equal(LC$rhs_list,list(Nmeas))
  
})
