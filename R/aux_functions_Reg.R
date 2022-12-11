
#############################################
#############################################
############         Reg         ############
#############################################
#############################################



ellipseCoeff <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
{
  A <- (t - k + 1) * eval_mean(cumSX, k, t)
  B <- (t - k + 1) * eval_mean(cumX, k, t)
  C <- (t - k + 1)
  D <- -(t - k + 1) * eval_mean(cumXY, k, t)
  E <- -(t - k + 1) * eval_mean(cumY, k, t)
  G <- (t - k + 1) * eval_mean(cumSY, k, t) + costQ[shift(k-1)] + beta

  return(c(A,B,C,D,E,G))
}

#############################################

evalReg_segment_kt <- function(cumX, cumY, cumXY, cumSX, cumSY, k, t)
{
  if(t == k){return(0)} #segment of size 1 => cost = +Inf (no such segment in regression)
  t1 <- (eval_mean(cumXY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  t2 <- (eval_mean(cumSX, k, t)*eval_mean(cumY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumXY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t) * t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 - eval_mean(cumSY, k, t))
  return(eval)
}

#############################################

evalReg_q_min <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta) ### minimum of q_{t}^{k}, data y_{k} to y_{t}
{
  if(t == k){return(Inf)}
  t1 <- (eval_mean(cumXY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  t2 <- (eval_mean(cumSX, k, t)*eval_mean(cumY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumXY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t)*t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 - eval_mean(cumSY, k, t)) + costQ[shift(k-1)] + beta
  return(eval)
}


evalReg_q_min2 <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta) ### minimum of q_{t}^{k}, data y_{k} to y_{t}
{
  if(t == k){return(costQ[shift(k-1)] + beta)} #### ATTENTION
  t1 <- (eval_mean(cumXY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  t2 <- (eval_mean(cumSX, k, t)*eval_mean(cumY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumXY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t)*t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 - eval_mean(cumSY, k, t)) + costQ[shift(k-1)] + beta
  return(eval)
}

#############################################

evalReg_q <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1, t2) ### value of q_{t}^{k}, data y_{k} to y_{t} at point (t1,t2)
{
  if(t == k){return(Inf)}
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t)*t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 - eval_mean(cumSY, k, t)) + costQ[shift(k-1)] + beta
  return(eval)
}


#############################################
#############################################
#############################################


evalReg_q_1_min <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t, beta)
{
  ### ### GET COEFF ### ###
  ### ### GET COEFF ### ###
  coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
  coeffk <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
  coeff1C <- coeffj - coeffk


  centerk <- ellipseCenter(coeffk)

  #print("WHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHY")
  #print(centerk)
  #print(ellipseEval(coeff1C, centerk$x, centerk$y))
  #print(ellipseEval(coeff1C, centerk$x, centerk$y) > 0)
  if(ellipseEval(coeff1C, centerk$x, centerk$y) > 0){return(-Inf)}
  #if(k == t){return(Inf)}
  #print("WHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHYWHY")


  ######################################################################
  #ell <- points.on.ellipse(coeff1C, 100)
  #minx <- min(c(ell[,1]))
  #maxx <- max(c(ell[,1]))
  #miny <- min(c(ell[,2]))
  #maxy <- max(c(ell[,2]))
  #plot(ell, asp = 1,  xlim = c(minx,maxx), ylim = c(miny,maxy))

  #print(c(j,k,t))
  #print(coeffk)
  #coeffk2 <- coeffk
  #A <- coeffk2[1]
  #B <- coeffk2[2]
  #C <- coeffk2[3]
  #print("uuuuuuuuuuuuuu")
  #print(A*C - B^2) # > 0
  #centerk <- ellipseCenter(coeffk2)
  #print(centerk)
  #GG <- ellipseEval(coeffk2, centerk$x, centerk$y)
  #coeffk2[6] <- coeffk2[6] - GG - 5
  #if(isAnEllipse(coeffk2))
  #{
  #  print(coeffk2)
  #  par(new = TRUE)
  #  plot(points.on.ellipse(coeffk2, 100), col = 2, asp = 1,  xlim = c(minx,maxx), ylim = c(miny,maxy))
  #}
  ######################################################################


  ### ### ROTATION MATRIX ### ###
  ### ### ROTATION MATRIX ### ###
  temp <- ellipseRotation(coeffk)
  angle <- temp$angle #rotation MATRIX


  ### ### ROTATION ### ###
  ### ### ROTATION ### ###
  coeffk <- temp$coeff
  temp <- ellipseRotation(coeff1C, angle)
  coeff1C <- temp$coeff


  ######################################################################
  #ell <- points.on.ellipse(coeff1C, 100)
  #minx <- min(c(ell[,1]))
  #maxx <- max(c(ell[,1]))
  #miny <- min(c(ell[,2]))
  #maxy <- max(c(ell[,2]))
  #plot(ell, asp = 1,  xlim = c(minx,maxx), ylim = c(miny,maxy))

  #print(c(j,k,t))
  #print(coeffk)
  #coeffk2 <- coeffk
  #A <- coeffk2[1]
  #B <- coeffk2[2]
  #C <- coeffk2[3]
  #print("uuuuuuuuuuuuuu")
  #print(A*C - B^2) # > 0
  #centerk <- ellipseCenter(coeffk2)
  #print(centerk)
  #GG <- ellipseEval(coeffk2, centerk$x, centerk$y)
  #coeffk2[6] <- coeffk2[6] - GG - 5
  #if(isAnEllipse(coeffk2))
  #{
  #  print(coeffk2)
  #  par(new = TRUE)
  #  plot(points.on.ellipse(coeffk2, 100), col = 2, asp = 1,  xlim = c(minx,maxx), ylim = c(miny,maxy))
  #}
  ######################################################################



  ### ### ALL ANGLES ### ###
  ### ### ALL ANGLES ### ###
  centerk <- ellipseCenter(coeffk)
  angle1 <- ellipseAngleFromConstant(coeff1C, centerk$x, "x")
  angle2 <- ellipseAngleFromConstant(coeff1C, centerk$y, "y")
  angles <- sort(c(unlist(angle1), unlist(angle2), c(0, pi/2, pi, 3*pi/2)))

  #print(angle1)
  #print(angle2)

  ### ### MIN SEARCH between two consecutive angles ### ###
  ### ### MIN SEARCH between two consecutive angles ### ###

  fk <- function(x) evalReg_q_onConstraint(coeffk, coeff1C, x)

  #print("anglesanglesanglesanglesanglesanglesanglesangles")
  #print(angles)
  #print(ellipseEval(coeff1C, centerk$x, centerk$y))
  #print(ellipseEval(coeff1C, centerk$x, centerk$y) > 0)

  l1k <- golden_Search(fk, angles[1], angles[2])
  l2k <- golden_Search(fk, angles[2], angles[3])
  l3k <- golden_Search(fk, angles[3], angles[4])
  l4k <- golden_Search(fk, angles[4], angles[5])
  l5k <- golden_Search(fk, angles[5], angles[6])
  l6k <- golden_Search(fk, angles[6], angles[7])
  l7k <- golden_Search(fk, angles[7], angles[8])
  l8k <- golden_Search(fk, angles[8], 2*pi)

  M <- c(l1k$m, l2k$m, l3k$m, l4k$m, l5k$m, l6k$m, l7k$m, l8k$m)
  return(min(M))
}

################################################
################################################

test <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t, beta)
{

  fk <- function(x) evalReg_q_onConstraint(coeffk, coeff1C, x)
  fkPrime <- function(x) evalReg_qPrime_onConstraint(coeffk, coeff1C, x)


  Newton_Raphson(fk, fkPrime, angles[1], angles[2])


}








