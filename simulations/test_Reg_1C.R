

data <- dataGenerator_Reg(c(100, 200), A = c(2,1), B = c(0,1), sdNoise = 1)

plot(data, col = c(rep(1,10),rep(2,10)))
beta <- 1*log(nrow(data))

(res1 <- OP_Reg(data, beta))
(res2 <- OP_Reg_PELT(data, beta))
(res3 <- OP_Reg_1C(data, beta))

all(res2$changepoints == res3$changepoints)


plot(res2$nb, type = 'b', col = 1, ylim = c(0, max(c(res3$nb,res2$nb))))
par(new = TRUE)
plot(res3$nb, type = 'b', col = 2, ylim = c(0, max(c(res3$nb,res2$nb))))


res1$changepoints
res2$changepoints
res3$changepoints

globalCost_Reg(data, res1$changepoints, beta)
globalCost_Reg(data, res2$changepoints, beta)
globalCost_Reg(data, res3$changepoints, beta)

res3$nb - res2$nb

res2$nb
res3$nb


