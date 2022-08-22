

n <- 1000
myBeta <- 4*log(n)
v12 <- 1
data <- dataGenerator2D(chpts <- 1:n,
                        means1 = rnorm(1000, sd = 1),
                        means2 = rnorm(1000, sd = 1),
                        sdNoise1 = v12,
                        sdNoise2 = v12)
plot(data$y2)
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



n <- 20000
myBeta <- 4*log(n)
v12 <- 1
data <- dataGenerator2D(chpts = c(n/2,n),
                        means1 = c(0,2),
                        means2 = c(0,1),
                        sdNoise1 = 1,
                        sdNoise2 = 1)
plot(data$y2)

system.time(R1 <- OP_2D_PELT(data = data, beta = myBeta))
system.time(R2 <- OP_2D_1C(data = data, beta = myBeta))
all(R1$changepoints == R2$changepoints)



which(R2$nb-R1$nb>0)

ymax <- log(max(R1$nb, R2$nb))
plot(log(R1$nb), type = 'l', ylim = c(0,ymax))
par(new = TRUE)
plot(log(R2$nb), type = 'l', col = 2, ylim = c(0,ymax))
R0 <- OP_2D(data,beta =myBeta)

R0$changepoints
R1$changepoints
R2$changepoints

globalCost(data, R1$changepoints, myBeta)
globalCost(data, R0$changepoints, myBeta)
globalCost(data, R2$changepoints, myBeta)

plot(R2$nb - R1$nb, type = 'l', col = 3)

all(R2$nrows/R2$nb == (R2$nb+1)/2)


