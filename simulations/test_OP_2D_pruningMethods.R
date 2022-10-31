
n <- 100
myBeta <- 1*log(n)
v12 <- 1
data <- dataGenerator2D(chpts = c(n),
                        means1 = c(1),
                        means2 = c(2),
                        sdNoise1 = 1,
                        sdNoise2 = 1)

plot(data$y2)
a <- system.time(S1 <- OP_2D_1C(data = data, beta = myBeta, testMode = 0))
b <- system.time(S2 <- OP_2D_1C(data = data, beta = myBeta, testMode = 1))
c <- system.time(S3 <- OP_2D_2C(data = data, beta = myBeta, testMode = 0))
d <- system.time(S4 <- OP_2D_2C(data = data, beta = myBeta, testMode = 1))
e <- system.time(S5 <- OP_2D_2C(data = data, beta = myBeta, testMode = 2))
f <- system.time(S6 <- OP_2D_AC(data = data, beta = myBeta))


all(S1$changepoints == S2$changepoints)
all(S1$changepoints == S3$changepoints)
all(S1$changepoints == S4$changepoints)
all(S1$changepoints == S5$changepoints)
all(S1$changepoints == S6$changepoints)

S1$changepoints
S6$changepoints

a[[1]]
b[[1]]
c[[1]]
d[[1]]
e[[1]]
f[[1]]

#####

n <- 500
myBeta <- 4*log(n)
v12 <- 1
data <- dataGenerator2D(chpts = c(n),
                        means1 = c(1),
                        means2 = c(2),
                        sdNoise1 = 1,
                        sdNoise2 = 1)


a <- system.time(S1 <- OP_2D_PELT(data = data, beta = myBeta))
b <- system.time(S2 <- OP_2D_1C(data = data, beta = myBeta, testMode = 0))
f <- system.time(S6 <- OP_2D_AC(data = data, beta = myBeta))

aa
a
f



S1
S6

S1$nb
S2$nb
S6$nb

S1$nrows
S2$nrows
S6$nrows

S1$changepoints
S2$changepoints
S6$changepoints


S1$nrows
S6$nrows

