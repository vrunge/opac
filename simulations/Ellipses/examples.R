

A <- 2
B <- 0.9
C <- 1
D <- -0.5
E <- 0.5
G <- -1.9

coeff <- c(A,B,C,D,E,G)
E1C <- my.ellipse(coeff)
cen <- ellipseCenter(coeff)
base::plot(E1C, col = 1, type = 'l', asp = 1, ylim = c(-2,2), xlim = c(-2,2))

E2 <- my.ellipse(c(1,0,1,0,0,-1))
cen <- ellipseCenter(c(1,0,1,0,0,-1))
par(new = TRUE)
base::plot(E2, col = 2, type = 'l', asp = 1, ylim = c(-2,2), xlim = c(-2,2))
points(cen$x, cen$y, pch = 12)

a <- Re((-E-2*D*1i)/(2*B-(C-A)*1i))
b <- Im((-E-2*D*1i)/(2*B-(C-A)*1i))
a
b
segments(cen$x, cen$y, 7*a, 7*b)

acos(a/b)

##############################
##############################
##############################

segments(cen$x, cen$y, a, -b)

theta <- seq(0,pi/2, length.out = 1000)
fun <- function(x)
{
  (C-A)*sin(2*x) + 2*B*cos(2*x) - 2*D*sin(x) +2*E*cos(x)
}

distance <- function(coeff, theta)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]
  u <- D*cos(theta) + E*sin(theta)
  v <- A*cos(theta)^2 + 2*B*cos(theta)*sin(theta) + C*sin(theta)^2
  return((-u +sqrt(u^2 - G*v))/v)
}

plot(theta,distance(coeff, theta), type = 'l')

plot(theta,fun(theta))
abline(v = acos(a/b))
