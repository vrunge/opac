

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



