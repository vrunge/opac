

#############################################

evalReg_q_1_min_Approx <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t, beta)
{
  ### ### GET COEFF ### ###
  coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
  coeffk <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
  coeff1C <- coeffj - coeffk

  centerj <- ellipseCenter(coeffj)
  centerk <- ellipseCenter(coeffk)
  center1C <- ellipseCenter(coeff1C)

  t1jt <- centerj$x
  t2jt <- centerj$y
  t1kt <- centerk$x
  t2kt <- centerk$y
  t1jk <- center1C$x
  t2jk <- center1C$y

  Dj <- sqrt((t1jt - t1jk)^2 + (t2jt - t2jk)^2)
  Dk <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)

  axes_1C <- ellipseAxes(coeff1C)
  R_small <- axes_1C$b

  A <- coeffk[1]
  B <- coeffk[2]
  C <- coeffk[3]
  Rk <- (1/2)*(A+C - sqrt((A-C)^2 + 4*B^2))
  A <- coeffj[1]
  B <- coeffj[2]
  C <- coeffj[3]
  Rj <- (1/2)*(A+C - sqrt((A-C)^2 + 4*B^2))

  minj <- evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
  mink <- evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)

  mj <- minj + Rj*(Dj - R_small)^2
  mk <- mink + Rk*(Dk - R_small)^2

  if(is.na(mj) | is.na(mk)){return(Inf)}
  res <- min(mj, mk)
  res <- trunc(res*1e13)*1e-13
  return(res)
}



