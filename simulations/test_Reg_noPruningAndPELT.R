


data <- dataGenerator_Reg(c(150, 300), A = c(2,1.5), B = c(0,1), sdNoise = 1)

plot(data, col = c(rep(1,150),rep(2,150)))

(res2 <- OP_Reg(data))
(res4 <- OP_Reg_PELT(data))

all(res2$changepoints == res4$changepoints)

plot(res4$nb)
res4

##############################

data <- dataGenerator_Reg(c(80, 150, 300), A = c(2,1.5,1), B = c(0,-1,1), sdNoise = 1)

plot(data)
plot(data, col = c(rep(1,80),rep(2,70),rep(3,150)))

(res2 <- OP_Reg(data))
(res4 <- OP_Reg_PELT(data))

all(res2$changepoints == res4$changepoints)

plot(res4$nb)
res4



##############################
###### OP_Reg + OP_Reg_PELT
##############################
data <- dataGenerator_Reg(chpts = c(50,100,150), A = c(-1,1,-1), B = c(-1,1,-1),
                          sdX = 5, sdNoise = 1)

Reg1 <- OP_Reg(data)
Reg2 <- OP_Reg_PELT(data)


### CHANGEPOINT

all(Reg1$changepoints == Reg2$changepoints)
all(Reg1$changepoints == Reg2$changepoints)

Reg1$changepoints

globalCost_Reg(data, Reg1$changepoints, 0)
globalCost_Reg(data, Reg2$changepoints, 0)




