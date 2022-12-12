
a = b = c = 0
i <- 1
while((a==b) & (b ==c))
{
  print(i)
  i <- i+1
  data <- dataGenerator_2D(chpts = 100)
  be <- runif(1,min = 1, max = 10)
  (res0 <- OP_2D(data, beta = be))
  (res1 <- OP_2D_PELT(data, beta = be))
  (res2 <- OP_2D_1C(data, beta = be))
  (res3 <- OP_2D_2C(data, beta = be))
  (res4 <- OP_2D_AC(data))
  res0$changepoints
  res1$changepoints
  res2$changepoints
  res3$changepoints
  res4$changepoints
  globalCost_2D(data,res0$changepoints, be)
  globalCost_2D(data,res1$changepoints, be)
  globalCost_2D(data,res2$changepoints, be)
  globalCost_2D(data,res3$changepoints, be)
  globalCost_2D(data,res4$changepoints, be)
  print(sum(res4$nb > res3$nb))
  res4$nb
  res3$nb
  a <- globalCost_2D(data,res2$changepoints, be)
  b <- globalCost_2D(data,res3$changepoints, be)
  c <- globalCost_2D(data,res4$changepoints, be)
}

a
b
c


###
#
#
#
#
#
#
#


dat <- dataGenerator_2D(chpts = 20)
b <- 2
(res1 <- OP_2D_PELT(dat, beta = b))
(res2 <- OP_2D_1C(dat, beta = b))
(res3 <- OP_2D_2C(dat, beta = b))
(res4 <- OP_2D_AC(dat, beta = b))
res1$changepoints
res2$changepoints
res3$changepoints
res4$changepoints




dat <- dataGenerator_2D(chpts = 5)
b <- 0
dat[1,] = dat[2,]
dat[3,] = dat[2,]
(res1 <- OP_2D_PELT(dat, beta = b))
(res2 <- OP_2D_1C(dat, beta = b))
(res3 <- OP_2D_2C(dat, beta = b))
(res4 <- OP_2D_AC(dat, beta = b))
res1$changepoints
res2$changepoints
res3$changepoints
res4$changepoints



