

##############################
##############################
##############################


n <- 100
myBeta <- 4*log(n)
v12 <- 1
data <- dataGenerator2D(chpts <- 1:n,
                        means1 = rnorm(n, sd = 1),
                        means2 = rnorm(n, sd = 1),
                        sdNoise1 = v12,
                        sdNoise2 = v12)
plot(data$y2)
R1 <- OP_2D(data,beta =myBeta)
R1$changepoints
R2 <- OP_2D_PELT(data,beta =myBeta)
R2$changepoints
all(R1$changepoints == R2$changepoints)
R1
R2

### plot ###
R2$changepoints
plot(R2$nb, type = 'l')
abline(v = R2$changepoints, col = 2)


############################################################################
############################################################################


n <- 100
myBeta <- 4*log(n)
v12 <- 1
data <- dataGenerator2D(chpts = c(n),
                        means1 = c(1),
                        means2 = c(2),
                        sdNoise1 = 1,
                        sdNoise2 = 1)

plot(data$y2)
R1 <- OP_2D(data = data, beta = myBeta)
R2 <- OP_2D_PELT(data = data, beta = myBeta)
R3 <- OP_2D_1C(data = data, beta = myBeta)


### CHANGEPOINT

all(R1$changepoints == R2$changepoints)
all(R1$changepoints == R3$changepoints)

R1$changepoints
R2$changepoints
R3$changepoints



globalCost_2D(data, R1$changepoints, myBeta)
globalCost_2D(data, R2$changepoints, myBeta)
globalCost_2D(data, R3$changepoints, myBeta)


### PRUNING strengh

######################################################



ymax <- (max(R2$nb, R3$nb))
plot((R2$nb), type = 'l', ylim = c(0,ymax))
par(new = TRUE)
plot((R3$nb), type = 'l', col = 2, ylim = c(0,ymax))


###

1 - ((R3$nrows-R3$nb) / (R3$nb*(R3$nb-1)/2)) ### reduction of the cross terms in percent
R3$nb

sum(1 - (R3$nrows-R3$nb) / (R3$nb*(R3$nb-1)/2) > 0.4, na.rm = T)






############################################################################
############################################################################
############################################################################
############################################################################
############################################################################


R4 <- OP_2D_2C(data = data, beta = myBeta)
R5 <- OP_2D_AC(data = data, beta = myBeta)





### CHANGEPOINT
all(R3$changepoints == R4$changepoints)
all(R3$changepoints == R5$changepoints)

R3$changepoints

R3$nb
R4$nb
R3$nb - R4$nb
R4$nb - R5$nb



R4$nb / R3$nb
R4$nrows3
R4$nrows3 - R5$nrows3

R4$nrows3
R5$nrows3




