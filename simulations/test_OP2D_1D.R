


data <- dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4))

res1 <- OP_2D_PELT(data, beta = 1)
res2 <- OP_2D_1C(data, beta = 1)
res1$changepoints
res2$changepoints
res0 <- OP_2D(data, beta = 6)
res1 <- OP_2D_PELT(data, beta = 6)
res2 <- OP_2D_1C(data, beta = 6)

res0$changepoints
res1$changepoints
res2$changepoints


res1 <- OP_2D_PELT(data, beta = 10)
res2 <- OP_2D_1C(data, beta = 10)
res1$changepoints
res2$changepoints
res1 <- OP_2D_PELT(data, beta = 1600)
res2 <- OP_2D_1C(data, beta = 1600)
res1$changepoints
res2$changepoints

for(i in 1:100)
{
data <- dataGenerator_2D(chpts = 500, means1 = 0, means2 = 0)

res1 <- OP_2D_PELT(data, 3)
res2 <- OP_2D_1C(data, 3)

print(globalCost_2D(data, res1$changepoints, 0) == globalCost_2D(data, res2$changepoints, 0))

res1 <- OP_2D_PELT(data, 6)
res2 <- OP_2D_1C(data, 6)

print(globalCost_2D(data, res1$changepoints, 0) == globalCost_2D(data, res2$changepoints, 0))


res1 <- OP_2D_PELT(data, 10)
res2 <- OP_2D_1C(data, 10)

print(globalCost_2D(data, res1$changepoints, 0) == globalCost_2D(data, res2$changepoints, 0))

res1 <- OP_2D_PELT(data, 20)
res2 <- OP_2D_1C(data, 20)

print(globalCost_2D(data, res1$changepoints, 0) == globalCost_2D(data, res2$changepoints, 0))


}










res1 <- OP_2D_PELT(data, beta = 10)
res2 <- OP_2D_1C(data, beta = 10)
res1$changepoints
res2$changepoints
res1$nb - res2$nb  #### this is bad news



res1 <- OP_2D_PELT(data, beta = 1600)
res2 <- OP_2D_1C(data, beta = 1600)
res1$changepoints
res2$changepoints





data <- dataGenerator_2D(chpts = 5, means1 = 0, means2 = 0)
data[1,] = data[2,]
data[3,] = data[2,]
res0 <- OP_2D(data, beta = 0)
res1 <- OP_2D_PELT(data, beta = 0)
res2 <- OP_2D_1C(data, beta = 0)
res0$changepoints
res1$changepoints
res2$changepoints

res1 <- OP_2D_PELT(data, beta = 6)
res2 <- OP_2D_1C(data, beta = 6)
res1$changepoints
res2$changepoints
res1 <- OP_2D_PELT(data, beta = 10)
res2 <- OP_2D_1C(data, beta = 10)
res1$changepoints
res2$changepoints
res1 <- OP_2D_PELT(data, beta = 1600)
res2 <- OP_2D_1C(data, beta = 1600)
res1$changepoints
res2$changepoints

###

data <- dataGenerator_2D(chpts = 1:500, means1 = rnorm(500,sd=10), means2 = rnorm(500,sd =10))

res1 <- OP_2D_PELT(data)
res2 <- OP_2D_1C(data)
res1$changepoints
res2$changepoints
res1$nb - res2$nb

res1 <- OP_2D_PELT(data, beta = 10)
res2 <- OP_2D_1C(data, beta = 10)
res1$changepoints
res2$changepoints
res1$nb - res2$nb



