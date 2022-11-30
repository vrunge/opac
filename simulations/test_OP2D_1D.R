


data <- dataGenerator_2D(chpts = 10*c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4))


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


data <- dataGenerator_2D(chpts = 500, means1 = 0, means2 = 0)

res1 <- OP_2D_PELT(data)
res2 <- OP_2D_1C(data)
res1$changepoints
res2$changepoints
res1 <- OP_2D_PELT(data, beta = 6)
res2 <- OP_2D_1C(data, beta = 6)
res1$changepoints
res2$changepoints
res1 <- OP_2D(data, beta = 10)
res2 <- OP_2D_1C(data, beta = 10)
res1$changepoints
res2$changepoints
res1 <- OP_2D_PELT(data, beta = 1600)
res2 <- OP_2D_1C(data, beta = 1600)
res1$changepoints
res2$changepoints





data <- dataGenerator_2D(chpts = 5, means1 = 0, means2 = 0)
#data[1,] = data[2,]
#data[3,] = data[2,]
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




