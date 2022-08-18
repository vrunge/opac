

n <- 1000
myBeta <- 10*log(n)
v12 <- 1
data <- dataGenerator2D(chpts <- 1:n,
                        means1 = rnorm(1000, sd = 1),
                        means2 = rnorm(1000, sd = 1),
                        sdNoise1 = v12,
                        sdNoise2 = v12)
R1 <- OP_2D(data,beta =myBeta)
R1$changepoints
R2 <- OP_2D_PELT(data,beta =myBeta)
R2$changepoints
all(R1$changepoints == R2$changepoints)


### plot ###
R2$changepoints
plot(R2$nb, type = 'l')
abline(v = R2$changepoints, col = 2)


############################################################################



n <- 1000
myBeta <- 1*log(n)
v12 <- 1
data <- dataGenerator2D(chpts = c(n/2,n),
                        means1 = c(0,2),
                        means2 = c(0,1),
                        sdNoise1 = 1,
                        sdNoise2 = 1)






R1 <- OP_2D_PELT(data = data, beta = myBeta)
R2 <- OP_2D_1C(data = data, beta = myBeta)
all(R1$changepoints == R2$changepoints)
R2$nb
R2$nrows
R1$changepoints
R2$changepoints


R1$nb
R1$costQ
R2$costQ

which(R2$nb[2:n]-R1$nb[1:(n-1)]>0)

ymax <- max(R1$nb, R2$nb)
plot(R1$nb[1:(n-1)], type = 'l', ylim = c(0,ymax))
par(new = TRUE)
plot(R2$nb[2:n], type = 'l', col = 2, ylim = c(0,ymax))
R0 <- OP_2D(data,beta =myBeta)

R0$changepoints
R1$changepoints
R2$changepoints

globalCost(data, R1$changepoints, myBeta)
globalCost(data, R0$changepoints, myBeta)
globalCost(data, R2$changepoints, myBeta)


