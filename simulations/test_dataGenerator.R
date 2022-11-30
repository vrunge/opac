

###########################
##### dataGenerator2D #####
###########################

myData <- dataGenerator2D(chpts = c(30,100,120), means1 = c(0,5,0), means2 = c(7,1,-4))


plot(myData$y1)
plot(myData$y2)


myData <- dataGenerator2D(chpts = c(300,1000,1200), means1 = c(1,5,10), means2 = c(7,1,4), type = "poisson")

plot(myData$y1)
plot(myData$y2)



###################################
##### dataGeneratorRegression #####
###################################

myData <- dataGeneratorRegression(chpts = c(40,90),
                                  A = c(2,-1),
                                  B = c(-10,2),
                                  meansX = c(10,2),
                                  sdX = c(1,10),
                                  sdNoise = 1)

plot(myData$x, myData$y, col = c(rep(1,40), rep(2,90)))


##############################
##### dataGeneratorSlope #####
##############################


myData <- dataGeneratorSlope(chpts = c(1,30,100,150), kinks = c(-2,pi,0.001,2.1), varNoise = 0)
myData
plot(myData)



myData <- dataGeneratorSlope(chpts = c(1,400,990,1500), kinks = c(0,4,-2,0), varNoise = 1)
myData
plot(myData)



myData <- dataGeneratorSlope(chpts = c(1,400,990,1500), kinks = c(1,10,5,15), varNoise = 1, type = "poisson")
myData
plot(myData)




