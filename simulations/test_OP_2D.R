

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


####



n <- 1000
myBeta <- 20*log(n)
v12 <- 1
data <- dataGenerator2D(chpts = n,
                        means1 = 0,
                        means2 = 0,
                        sdNoise1 = 2,
                        sdNoise2 = 1)



R1 <- OP_2D_PELT(data = data, beta = myBeta)
all(R1$changepoints == R2$changepoints)

R2 <- OP_2D_1C(data = data, beta = myBeta)
R2$nb
R2$nrows

ymax <- max(R1$nb, R2$nb)
plot(R1$nb, type = 'l', ylim = c(0,ymax))
par(new = TRUE)
plot(R2$nb, type = 'l', col = 2, ylim = c(0,ymax))


plot(R2$nb-R1$nb, type = 'l', col = 2)




