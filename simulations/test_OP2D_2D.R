

data <- dataGenerator_2D(chpts = 200)

res1 <- OP_2D_PELT(data)
res2 <- OP_2D_1C(data)
res3 <- OP_2D_2C(data)
res1$changepoints
res2$changepoints
res3$changepoints

res2
res3

res2$nb - res3$nb
