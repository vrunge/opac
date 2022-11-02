

#### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### inscribed and circumscribed circles from Ellipses #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### ####

#coeff Ax^2 + Bxy + Cy^2 + Dx + Ey + F
coeff <- c(10,-5,2,4,2,-6)

E1 <- my.ellipse(coeff)
Ein <- my.ellipse(new.coeff.for.circle(coeff, "in"))
Eout <- my.ellipse(new.coeff.for.circle(coeff, "out"))
minx <- min(E1[,1])-1.5
maxx <- max(E1[,1])+1.5
miny <- min(E1[,2])-1.5
maxy <- max(E1[,2])+1.5

base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1)
par(new = TRUE)
base::plot(Ein, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1)
par(new = TRUE)
base::plot(Eout, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 3, type = 'l', asp = 1)



#############################################################################
######### from data x and y

n <- 10
eps <- 1
a <- 2
b <- 3
x <- rnorm(n,0.3)
y <- a*x + b + rnorm(n,sd = eps)


res <- my.coeff.ellipse.from.data(x, y, 100)
res


E1 <- my.ellipse(res$coeffEllipse)
Ein <- my.ellipse(res$coeffC1)
Eout <- my.ellipse(res$coeffC2)

minx <- min(E1[,1])-1.5
maxx <- max(E1[,1])+1.5
miny <- min(E1[,2])-1.5
maxy <- max(E1[,2])+1.5
base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1)
par(new = TRUE)
base::plot(Ein, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1)
par(new = TRUE)
base::plot(Eout, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 3, type = 'l', asp = 1)

res





#### #### #### #### #### #### #### #### #### ####
#### #### Problem B with one constraint #### ####
#### #### #### #### #### #### #### #### #### ####


n <- 100
eps <- 1
a <- 2
b <- 3
x <- rnorm(n,0)
y <- a*x + b + rnorm(n,sd = eps)

res1 <- my.coeff.ellipse.from.data(x, y, 1000)
x <- rnorm(n,0)
y <- a*x + b + rnorm(n,sd = eps)

res2 <- my.coeff.ellipse.from.data(x, y, 1000)


E1 <- my.ellipse(res1$coeffEllipse)
E2 <- my.ellipse(res2$coeffEllipse)

minx <- min(E1[,1])-1.5
maxx <- max(E1[,1])+1.5
miny <- min(E1[,2])-1.5
maxy <- max(E1[,2])+1.5
base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1)
par(new = TRUE)
base::plot(E2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1)



optim <- dualFunctionB_q_1(res1$coeffEllipse, res2$coeffEllipse, seq(-2,2,0.01))

plot(optim$t1,optim$t2, type = 'b')

plot(optim$m, type = 'l')
optim

###

E1 <- my.ellipse(c(1,-0.5,1,1,1,-5))
E2 <- my.ellipse(c(1,1,3,2,2,-5))
minx <- min(E1[,1])-1.5
maxx <- max(E1[,1])+1.5
miny <- min(E1[,2])-1.5
maxy <- max(E1[,2])+1.5
base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1)
par(new = TRUE)
base::plot(E2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1)



optim <- dualFunctionB_q_1(c(1,1,1,1,1,-5), c(1,-2,1,2,2,-5), seq(-2,2,0.01))

plot(optim$t1,optim$t2, type = 'b')

plot(optim$m, type = 'l')


index <- which.max(optim$m)
index

base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1)
par(new = TRUE)
base::plot(E2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1)
points(optim$t1[index],optim$t2[index])


###################################

c1 <- c(1,0,1,1,1,-5)
c2 <- c(1,0,3,2,2,-10)
optim <- dualFunctionB(c1, c2, seq(-5,5,0.01))

plot(optim$t1,optim$t2, type = 'b')

plot(optim$m, type = 'l')


index <- which.max(optim$m)
index

E1 <- my.ellipse(c1)
E2 <- my.ellipse(c2)
base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1)
par(new = TRUE)
base::plot(E2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1)
points(optim$t1[index],optim$t2[index])



#### # # #
optim <- primalFunctionB(c(1,0,1,1,1,-5), c(1,0,1,2,2,-7), seq(-5,5,0.01))
plot(optim$p, type = 'l')







