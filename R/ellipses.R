

isAnEllipse <- function(coeff)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  if(abs(G) == Inf){return(FALSE)}
  u <- A*C - B^2
  v <- B*D - A*E
  w <- v^2 + u*(D^2 - A*G)
  if((u < 10^(-14)) | (w < 10^(-14))){return(FALSE)}
  return(TRUE)
}

#############################################

ellipseCenter <- function(coeff)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]
  u <- A*C - B^2 # > 0
  #NULL (one of the minimal points)
  #translated paraboloid
  if(u < 10^(-14)){return(list(x = 0, y = -E/C))}
  if(u < 0){return(NULL)} #NULL
  x <- (B*E - C*D)/u
  y <- (B*D - A*E)/u
  return(list(x = x, y = y))
}

#############################################

ellipseAngle <- function(coeff)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  if((B == 0) && (A <= C)){return(0)}
  if((B == 0) && (A > C)){return(pi/2)}
  value <- (1/(2*B))*(C-A - sqrt((C-A)^2 + 4*B^2))
  return(atan(value))
}

#############################################

ellipseAxes <- function(coeff)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]
  u <- A*C - B^2
  va <- A+C + sqrt((A-C)^2 + 4*B^2)
  vb <- A+C - sqrt((A-C)^2 + 4*B^2)
  w <- 2*(A*E^2 + C*D^2 - 2*B*D*E + (B^2 - A*C)*G)
  a <- (1/(2*u))*sqrt(w)*sqrt(va)
  b <- (1/(2*u))*sqrt(w)*sqrt(vb)
  return(list(a = a, b = b))
}

#############################################

ellipseRotation <- function(coeff, angle = NULL)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]
  M <- matrix(c(A,B,B,C),2,2)

  if(is.null(angle))
  {
    angle <- atan((1/(2*B))*(C - A - sqrt((A-C)^2 + 4*B^2)))
  }
  rot <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2,2)
  newM <- t(rot)%*%M%*%rot
  newDE <- (t(rot))%*%c(D,E)
  newA <- newM[1,1]
  newB <- newM[1,2]
  newC <- newM[2,2]
  newD <- newDE[1]
  newE <- newDE[2]
  newG <- G
  return(list(coeff = c(newA, newB, newC, newD, newE, newG), angle = angle))
}

#############################################

ellipseEval <- function(coeff, t1, t2)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]
  eval <- A*t1^2 + 2*B*t1*t2 + C*t2^2 + 2*D*t1 + 2*E*t2 + G
  return(eval)
}

#############################################

ellipseEvalPrime <- function(coeff, t1, t2, d1, d2)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  eval <- A*d1*t1 + B*(d1*t2 + d2*t1) + C*d2*t2 + 2*D*d1 + 2*E*d2
  return(2*eval)
}

#############################################

ellipseCoordinatesFromAngle <- function(coeff, angle)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  center <- ellipseCenter(coeff)
  theta <- ellipseAngle(coeff)
  axes <- ellipseAxes(coeff)

  x <- axes$a*cos(theta)*cos(angle) - axes$b*sin(theta)*sin(angle) + center$x
  y <- axes$a*sin(theta)*cos(angle) + axes$b*cos(theta)*sin(angle) + center$y

  return(list(x = x, y = y))
}

#############################################


ellipseDerivativesFromAngle <- function(coeff, angle)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  center <- ellipseCenter(coeff)
  theta <- ellipseAngle(coeff)
  axes <- ellipseAxes(coeff)

  x <- -axes$a*cos(theta)*sin(angle) - axes$b*sin(theta)*cos(angle)
  y <- -axes$a*sin(theta)*sin(angle) + axes$b*cos(theta)*cos(angle)

  return(list(x = x, y = y))
}


#############################################

evalReg_q_onConstraint <- function(coeff_eval, coeff_constraint, theta)
{
  pointToEval <- ellipseCoordinatesFromAngle(coeff_constraint, theta)
  return(ellipseEval(coeff_eval, pointToEval$x, pointToEval$y))
}

#############################################

evalReg_qPrime_onConstraint <- function(coeff_eval, coeff_constraint, theta)
{
  pointToEval <- ellipseCoordinatesFromAngle(coeff_constraint, theta)
  derivatives <- ellipseDerivativesFromAngle(coeff_constraint, theta)
  return(ellipseEvalPrime(coeff_eval, pointToEval$x, pointToEval$y, derivatives$x, derivatives$y))
}

#############################################

ellipseAngleFromConstant <- function(coeff, const, type = "x")
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  center <- ellipseCenter(coeff)
  theta <- ellipseAngle(coeff)
  axes <- ellipseAxes(coeff)

  if(type == "x")
  {
    xbarre1 <- const
    xbarre2 <- const
    d2 <- (B*const + E)^2 - C*(A*const^2 + 2*D*const + G)
    if(d2 < 0){return(NULL)}
    ybarre1 <- -(B*const + E)/C - sqrt(d2)/C
    ybarre2 <- -(B*const + E)/C + sqrt(d2)/C
  }
  if(type == "y")
  {
    ybarre1 <- const
    ybarre2 <- const
    d2 <- (B*const + D)^2 - A*(C*const^2 + 2*E*const + G)
    if(d2 < 0){return(NULL)}
    xbarre1 <- -(B*const + D)/A - sqrt(d2)/A
    xbarre2 <- -(B*const + D)/A + sqrt(d2)/A
  }
  cos_u1 <- (1/(axes$a))*((xbarre1 - center$x)*cos(theta) + (ybarre1 - center$y)*sin(theta))
  sin_u1 <- (1/(axes$b))*(-(xbarre1 - center$x)*sin(theta) + (ybarre1 - center$y)*cos(theta))

  cos_u2 <- (1/(axes$a))*((xbarre2 - center$x)*cos(theta) + (ybarre2 - center$y)*sin(theta))
  sin_u2 <- (1/(axes$b))*(-(xbarre2 - center$x)*sin(theta) + (ybarre2 - center$y)*cos(theta))

  u1 <- atan2(sin_u1, cos_u1)
  u2 <- atan2(sin_u2, cos_u2)

  u1 <- u1 + 2*pi*(u1 < 0)
  u2 <- u2 + 2*pi*(u2 < 0)
  return(list(u1 = u1, u2 = u2))
}



#############################################



#' points.on.ellipse
#'
#' @description Generating 2d points describing an ellipse from its coefficients in Cartesian form.
#' @param coeff coefficients (A,B,C,D,E,F) in equation Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + F = 0
#' @param n.points number of (x,y) points to generate
#' @return matrix of size "2*n.points x 2" for (x,y) bipoints on the 2d plane describing an ellipse
points.on.ellipse <- function(coeff, n.points = 1000)
{
  res <- matrix(4*n.points, 2*n.points, 2)
  s <- seq(0, 1, length.out = n.points)
  for(i in 1: n.points)
  {
    temp <- ellipse_point_s_type(coeff, s[i], "right")
    res[i,] <- unlist(temp)
  }
  s <- seq(1, 0, length.out = n.points)
  for(i in 1: n.points)
  {
    temp <- ellipse_point_s_type(coeff, s[i], "left")
    res[n.points + i,] <- unlist(temp)
  }

  ###
  ### y in ]yMin, yMax[ to have two roots in x (for ellipses)
  ###

  return(res)
}


ellipse_point_s_type <- function(coeff, s, type)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  ###
  ### y in ]yMin, yMax[ to have two roots in x (for ellipses)
  ###
  u <- A*C - B^2
  v <- B*D - A*E
  w <- v^2 + u*(D^2 - A*G)

  if((u < 0) | (w < 0)){return(NULL)}
  if((s < 0) | (s > 1)){return(NULL)}

  yMin <- (v-sqrt(w))/u
  yMax <- (v+sqrt(w))/u

  ###
  ### computing y value
  ###
  y <- yMin + s*(yMax - yMin)

  if(abs(s - 1) < 10^(-14)){y <- yMax}
  if(abs(s) < 10^(-14)){y <- yMin}

  ###
  ### computing x value
  ###
  z1 <- B*y+D
  z2 <- z1^2 - A*(C*y^2 + 2*E*y + G)
  if((z2 < 10^(-13)) | (abs(s - 1)<10^(-14)) | (abs(s + 1)<10^(-14))){x <- -z1/A}
  else
  {
    if(type == "left"){x <- (-z1 - sqrt(z2))/A}
    if(type == "right"){x <- (-z1 + sqrt(z2))/A}
  }

  return(list(x = x, y = y))
}
