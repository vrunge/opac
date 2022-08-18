

n <- 1000
myBeta <- 10*log(n)
v12 <- 1
data <- dataGenerator2D(chpts <- 1:n,
                        means1 = rnorm(1000, sd = 3),
                        means2 = rnorm(1000, sd = 2),
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
