


data <- dataGenerator_Reg(c(150, 300), A = c(2,1.5), B = c(0,1), sdNoise = 1)

plot(data, col = c(rep(1,150),rep(2,150)))

(res1 <- OP_Reg(data))
(res2 <- OP_Reg_PELT(data))
(res3 <- OP_Reg_1C(data))

all(res2$changepoints == res4$changepoints)

plot(res3$nb)
res2
res3

