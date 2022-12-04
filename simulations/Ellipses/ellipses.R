

#' my.ellipse
#'
#' @description Generating 2d points describing an ellipse from its coefficients in Cartesian form.
#' @param coeff coefficients (A,B,C,D,E,F) in equation Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
#' @param n.points number of (x,y) points to generate
#' @return matrix of size "n.points x 2" for (x,y) bipoints on the 2d plane describing an ellipse
#' @examples
#' my.ellipse(c(10,-5,2,4,2,-6), 1000)
my.ellipse <- function(coeff, n.points = 1000)
{
  a <- coeff[1]
  b <- coeff[3]
  c <- coeff[2]
  d <- coeff[4]
  e <- coeff[5]
  f <- coeff[6]

  ## solve for centre
  A <- matrix(c(a, c / 2, c / 2, b), 2L)
  B <- c(-d / 2, -e / 2)
  mu <- solve(A, B)

  ## generate points on circle
  r <- sqrt(a * mu[1] ^ 2 + b * mu[2] ^ 2 + c * mu[1] * mu[2] - f)
  theta <- seq(0, 2 * pi, length = n.points)
  v <- rbind(r * cos(theta), r * sin(theta))

  z <- backsolve(chol(A), v) + mu ## transform for points on ellipse

  return(t(z))
}


#' new.coeff.for.circle
#'
#' @description Building coefficients for inscribed and circumscribed circles to an ellipse whose coefficients are given in Cartesian form.
#' @param coeff coefficients (A,B,C,D,E,F) in equation Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
#' @param type "in" or "out"
#' @return Six coefficients for Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0 (A=C and B=0)
#' @examples
#' new.coeff.for.circle(c(10,-5,2,4,2,-6), "in")
new.coeff.for.circle <- function(coeff, type = "in")
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  det <- B^2 - 4*A*C
  u <- 2*(A*E^2 + C*D^2-B*D*E+(B^2-4*A*C)*G)
  v <- A+C+sqrt((A-C)^2+B^2)
  w <- A+C-sqrt((A-C)^2+B^2)
  x0 <- (2*C*D-B*E)
  y0 <- (2*A*E-B*D)
  A2 <- det^2
  B2 <- 0
  C2 <- det^2
  D2 <- -2*x0*det
  E2 <- -2*y0*det
  if(type == "in"){G2 <- x0^2+y0^2-u*w}
  if(type == "out"){G2 <- x0^2+y0^2-u*v}
  return(c(A2,B2,C2,D2,E2,G2))
}



#' my.coeff.ellipse.from.data
#'
#' @description Building coefficients for inscribed and circumscribed circles to an ellipse from real vector data (x,y)
#' @param x data vector dim 1
#' @param y data vector dim 2 (linked to x by simple linear regression)
#' @param level a value equal to m_{k_1} - m_{j_1}
#' @return A list of elements
#' @examples
#' new.coeff.for.circle(c(10,-5,2,4,2,-6), "in")
my.coeff.ellipse.from.data <- function(x, y, level)
{
  n <- length(x)
  Vx <- ((n-1)/n)*var(x)
  Vy <- ((n-1)/n)*var(y)
  corxy <- (mean(x*y) - mean(x)*mean(y))/(sqrt(Vx*Vy))

  A <- sum(x^2)
  B <- 2*sum(x)
  C <- n
  D <- -2*sum(x*y)
  E <- -2*sum(y)
  G <- sum(y^2) - level
  coeffEllipse <- c(A,B,C,D,E,G)

  u <- (mean(x^2)+1)/Vx
  d1 <- (1/2)*(u + sqrt(u^2-4/Vx))
  d2 <- (1/2)*(u - sqrt(u^2-4/Vx))
  D <- level/n - Vy*(1-corxy^2)

  lmxy <- lm(y~x)
  a <- lmxy$coefficients[2]
  b <- lmxy$coefficients[1]

  coeffC1 <- c(1,0,1,-2*a,-2*b,-D*d1 + a^2 + b^2)
  coeffC2 <- c(1,0,1,-2*a,-2*b,-D*d2 + a^2 + b^2)

  s <- 4*Vx/(A/n+1)^2
  e <- sqrt(1  - (1 - sqrt(1-s))/(1 + sqrt(1-s)))
  ratio <- (1-mean(x^2))/mean(x)
  xx <- ratio- sqrt(ratio^2+4)

  return(list(coeffEllipse = coeffEllipse,
              coeffC1 = coeffC1,
              coeffC2 = coeffC2,
              ratioInExcentricity = s,
              excentricity = e,
              ratio = ratio,
              xx = xx,
              cosangle = 1/sqrt(1+xx^2),
              sinangle = xx/sqrt(1+xx^2),
              angle = 1/2*atan(B/(A-C))/pi))
}



#' dualFunctionB_q_1
#'
#' @description Building the dual function in problem B with one constraint
#' @param coeff1 coefficients (A,B,C,D,E,F) in objective function Ax^2 + Bxy + Cy^2 + Dx + Ey + F
#' @param coeff2 coefficients (A,B,C,D,E,F) in constraint Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0
#' @param x vector x for plot
#' @return the dual function to problem B with one constraint
#' @examples
#' dualFunctionB_q_1(c(10,-5,2,4,2,-6), c(10,-5,2,4,2,-6), seq(-5,5,0.01))
dualFunctionB_q_1 <- function(coeff1, coeff2, x)
{
  coeff1 <- coeff1 * c(1,2,1,2,2,1)
  coeff2 <- coeff2 * c(1,2,1,2,2,1)
  a <- coeff1[1] + x*coeff2[1]
  b <- coeff1[2] + x*coeff2[2]
  c <- coeff1[3] + x*coeff2[3]
  d <- coeff1[4] + x*coeff2[4]
  e <- coeff1[5] + x*coeff2[5]
  f <- coeff1[6] + x*coeff2[6]
  t1 <- (b*e-c*d)/(a*c-b^2)
  t2 <- (b*d-a*e)/(a*c-b^2)
  m <- (-c*d^2 - a*e^2 + 2*b*d*e)/(a*c-b^2) + f
  return(list(t1= t1, t2 = t2, m = m, det = a*c-b^2))
}





dualFunctionB <- function(coeff1, coeff2, x)
{
  coeff1 <- coeff1 * c(1,2,1,2,2,1)
  coeff2 <- coeff2 * c(1,2,1,2,2,1)
  a <- coeff1[1] + x*coeff2[1]
  b <- coeff1[2] + x*coeff2[2]
  c <- coeff1[3] + x*coeff2[3]
  d <- coeff1[4] + x*coeff2[4]
  e <- coeff1[5] + x*coeff2[5]
  f <- coeff1[6] + x*coeff2[6]
  t1 <- -d/a
  t2 <- -e/c
  m <- -d^2/a -e^2/c + f
  return(list(t1= t1, t2 = t2, m = m, det = a*c-b^2))
}





primalFunctionB <- function(coeff1, coeff2, x)
{
  coeff1 <- coeff1 * c(1,2,1,2,2,1)
  coeff2 <- coeff2 * c(1,2,1,2,2,1)
  a <- coeff1[1] + x*coeff2[1]
  b <- coeff1[2] + x*coeff2[2]
  c <- coeff1[3] + x*coeff2[3]
  d <- coeff1[4] + x*coeff2[4]
  e <- coeff1[5] + x*coeff2[5]
  f <- coeff1[6] + x*coeff2[6]
  t1 <- (b*e-c*d)/(a*c-b^2)
  t2 <- (b*d-a*e)/(a*c-b^2)
  t1 <- (b*e-c*d)
  t2 <- (b*d-a*e)

  res <- coeff2[1]*t1^2 + 2*coeff2[2]*t1*t2 + coeff2[3]*t2^2 + coeff2[4]*t1*(a*c-b^2) + coeff2[5]*t2*(a*c-b^2) + coeff2[6]*(a*c-b^2)^2
  return(list(p = res))
}





