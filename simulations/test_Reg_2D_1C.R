


data <- dataGenerator_Reg(c(150, 300), A = c(2,1.5), B = c(0,1), sdNoise = 1)

plot(data, col = c(rep(1,150),rep(2,150)))
beta <- 1
(res1 <- OP_Reg(data))
(res2 <- OP_Reg_PELT(data))
(res3 <- OP_Reg_1C(data))

all(res2$changepoints == res3$changepoints)

plot(res3$nb)
res1$changepoints
res2$changepoints
res3$changepoints

globalCost_Reg(data, res1$changepoints, beta)
globalCost_Reg(data, res2$changepoints, beta)
globalCost_Reg(data, res3$changepoints, beta)


