
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
  E <- (t - k + 1) * eval_mean(cumY, k, t)
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

#############################################

evalReg_q <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1, t2) ### minimum of q_{t}^{k}, data y_{k} to y_{t}
{
  if(t == k){return(Inf)}
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t)*t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 - eval_mean(cumSY, k, t)) + costQ[shift(k-1)] + beta
  return(eval)
}

#############################################


evalReg_q_1_min <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t, beta)
{
  coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
  coeffk <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
  coeff1C <- coeffj - coeffk

  center <- ellipseCenter(coeffk)

  f <- function(s, type)
  {
    branch <- ellipseBranch(coeff1C, s, type = type, ref0 = center$y)
    res <- evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, branch$x, branch$y)
    return(res)
  }

  # 2 explorations s in [0, 1] left and right
  # 2 explorations s in [-1, 0] left and right
  l1 <- golden_Search(f, 0, 1, "left")
  l2 <- golden_Search(f, -1, 0, "left")
  r1 <- golden_Search(f, 0, 1, "right")
  r2 <- golden_Search(f, -1, 0, "right")

  i <- which.min(c(l1$m, l2$m, r1$m, r2$m))
  m <- c(l1$m, l2$m, r1$m, r2$m)[i]
  s <- c(l1$s, l2$s, r1$s, r2$s)[i]
  p <- ellipseBranch(coeff1C, s, type = c("left", "left", "right", "right")[m], ref0 = center$y)

  return(list(p = c(p$x, p$y), m =  trunc(m*1e13)*1e-13))
}



##############################################
#####
##### Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + G #####
#####
##### s in [-1,1] with ref0 the return value y for which s = 0
#####
##############################################

ellipseBranch <- function(A, B, C, D, E, G, s, type, ref0)
{
  ###
  ### y in ]yMin, yMax[ to have two roots in x (for ellipses)
  ###
  u <- B^2 - A*C
  v <- B*D - A*E
  w <- v^2 - u*(D^2 - A*G)

  if((u > 0) | (w < 0)){return(NULL)}
  if((s < -1) | (s > 1)){return(NULL)}

  yMin <- (-v-sqrt(w))/u
  yMax <- (-v+sqrt(w))/u

  if(s >= 0){y <- ref0 + s*(yMax - ref0)}
  if(s < 0){y <- ref0 + s*(ref0 - yMin)}

  ###
  ### computing x value
  ###
  z1 <- B*y+D
  z2 <- z1^2 - A*(C*y^2 + 2*E*y + G)

  if(type == "left"){x <- (-z1 - sqrt(z2))/A}
  if(type == "right"){x <- (-z1 + sqrt(z2))/A}
  return(list(x = x, y = y))
}

#############################################

ellipseCenter <- function(A, B, C, D, E, G)
{
  u <- B^2 - A*C
  if(u > 0){return(NULL)}
  x <- (C*D-B*E)/u
  y <- (A*E-B*D)/u
  return(list(x = x, y = y))
}

#############################################

golden_Search <- function(f, a, b, type)
{
  phi <- 2/(sqrt(5) + 1)

  x1 <- b - phi*(b - a)
  x2 <- a + phi*(b - a)
  f1 <- f(x1, type)
  f2 <- f(x2, type)

  while (abs(b - a) > 10^(-14))
  {
    if (f2 > f1)
    {
      b <- x2
      x2 <- x1
      f2 <- f1
      x1 <- b - phi*(b - a)
      f1 <- f(x1, type)
    }
    else
    {
      a <- x1
      x1 <- x2
      f1 <- f2
      x2 <- a + phi*(b - a)
      f2 <- f(x2, type)
    }
  }
  return(list(s = (x1+x2)/2, m = f((x1+x2)/2)))
}
