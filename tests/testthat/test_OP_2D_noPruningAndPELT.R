
library(testthat)
library(opac)

test_that("2D change-point model. PELT is equal to OP", {
  data <- dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4))

  res1 <- OP_2D(data)
  res2 <- OP_2D_PELT(data)

  res2$changepoints
  expect_equal(res1$changepoints, res2$changepoints)
})
