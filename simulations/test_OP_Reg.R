

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
