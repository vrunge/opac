

##############################
###### OP_Reg + OP_Reg_PELT
##############################
data <- dataGeneratorRegression(chpts = c(50,100,150), A = c(-1,1,-1), B = c(-1,1,-1),
                                sdX = 5, sdNoise = 1)
Reg1 <- OP_Reg(data)
Reg2 <- OP_Reg_PELT(data)


### CHANGEPOINT

all(Reg1$changepoints == Reg2$changepoints)
all(Reg1$changepoints == Reg2$changepoints)

Reg1$changepoints

globalCost_Reg(data, Reg1$changepoints, myBeta)
globalCost_Reg(data, Reg2$changepoints, myBeta)


##############################
##############################
##############################


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
R1
R2

### plot ###
R2$changepoints
plot(R2$nb, type = 'l')

abline(v = R2$changepoints, col = 2)


############################################################################
############################################################################


n <- 100
myBeta <- 1*log(n)
v12 <- 1
data <- dataGenerator2D(chpts = c(n/2,n),
                        means1 = c(0,1),
                        means2 = c(0,2),
                        sdNoise1 = 1,
                        sdNoise2 = 1)

plot(data$y2)
R1 <- OP_2D(data = data, beta = myBeta)
R2 <- OP_2D_PELT(data = data, beta = myBeta)
R3 <- OP_2D_1C(data = data, beta = myBeta)


### CHANGEPOINT

all(R1$changepoints == R2$changepoints)
all(R1$changepoints == R2$changepoints)
all(R1$changepoints == R3$changepoints)

R1$changepoints
R3$changepoints



globalCost_2D(data, R1$changepoints, myBeta)
globalCost_2D(data, R2$changepoints, myBeta)
globalCost_2D(data, R3$changepoints, myBeta)


### PRUNING strengh

######################################################




R1$nb
(R2$nrows)/(R2$nb^2)/R2$nb




which(R2$nb-R1$nb>0)

ymax <- (max(R1$nb, R2$nb))
plot((R2$nb), type = 'l', ylim = c(0,ymax))
par(new = TRUE)
plot((R3$nb), type = 'l', col = 2, ylim = c(0,ymax))

R1$changepoints
R2$changepoints


plot(R3$nb - R2$nb, type = 'l', col = 3)

all(R2$nrows == R2$nb*(R2$nb+1)/2)

1 - (R3$nrows-R3$nb) / (R3$nb*(R3$nb-1)/2) ### reduction of the cross terms in percent
R3$nb

sum(1 - (R3$nrows-R3$nb) / (R3$nb*(R3$nb-1)/2) > 0.4, na.rm = T)

R2$nb
R3$nb

R2$changepoints
R3$changepoints

