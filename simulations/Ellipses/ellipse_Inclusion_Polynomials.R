

### inclusion => 3 racines réelles positives?
###################################
#coeff Ax^2 + Bxy + Cy^2 + Dx + Ey + F
c1 <- c(2,0,1,1,1,-10)
c2 <- c(2,0,3,2,2,-2)

polynom <- function(x, c1, c2)
{
  M <- matrix(0, 3, 3)
  M[upper.tri(M, diag=TRUE)] <- c1
  M1 <- (M + t(M))/2
  M <- matrix(0, 3, 3)
  M[upper.tri(M, diag=TRUE)] <- c2
  M2 <- (M + t(M))/2
  return(det(M1*x-M2))
}

x <- seq(0,3,length.out = 1000)
res <- sapply(x, function(x) polynom(x,c1 = c1,c2 = c2))

# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #




## PLOT
E1 <- my.ellipse(c1)
E2 <- my.ellipse(c2)
minx <- min(E1[,1])-1
maxx <- max(E1[,1])+1
miny <- min(E1[,2])-1
maxy <- max(E1[,2])+1
par(mfrow=c(1,2))
base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1, lwd = 2, xlab = "", ylab = "")
par(new = TRUE)
base::plot(E2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1, lwd = 2, xlab = "", ylab = "")
points(optim$t1[index],optim$t2[index])

plot(x,res, type = 'l', col = 4, lwd = 2, xlab = "", ylab = "")
abline(h = 0, v = 0,lwd = 2)




# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # #




## PLOT
c1 <- c(2,0,1,1,1,-10)
c2 <- c(2,0,3,2,2,-2)
par(mfrow=c(1,1))
E1 <- my.ellipse(c1)
E2 <- my.ellipse(c2)

minx <- min(E1[,1])-1
maxx <- max(E1[,1])+1
miny <- min(E1[,2])-1
maxy <- max(E1[,2])+1
base::plot(E1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 1, type = 'l', asp = 1, lwd = 2, xlab = "", ylab = "")
par(new = TRUE)
base::plot(E2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = 2, type = 'l', asp = 1, lwd = 2, xlab = "", ylab = "")

n <- 5
for(i in  1:n)
{
  u <- runif(1,0,4)
  c3 <- c1*u - c2
  print(c3)
  M <- matrix(0, 3, 3)
  M[upper.tri(M, diag=TRUE)] <- c3
  M3 <- (M + t(M))/2
  print(u)
  print(M3[1,1])
  print(det(M3[1:2,1:2]))
  print(det(M3))
  if(M3[1,1] > 0 & det(M3[1:2,1:2]) > 0 & det(M3) < 0)
  {
    E3 <-  my.ellipse(c1*u - c2)
    par(new = TRUE)
    base::plot(E3, xlim = c(minx,maxx), ylim = c(miny,maxy), col = i+2, type = 'l', asp = 1, lwd = 1, xlab = "", ylab = "")
  }
}



