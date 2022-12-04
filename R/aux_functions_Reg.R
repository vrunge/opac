
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
  if(t == k){return(0)}
  t1 <- (eval_mean(cumXY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  t2 <- (eval_mean(cumSX, k, t)*eval_mean(cumY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumXY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t)*t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 + eval_mean(cumSY, k, t))
  return(eval)
}

#############################################

evalReg_q_min <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta) ### minimum of q_{t}^{k}, data y_{k} to y_{t}
{
  if(t == k){return(Inf)}
  t1 <- (eval_mean(cumXY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  t2 <- (eval_mean(cumSX, k, t)*eval_mean(cumY, k, t) - eval_mean(cumX, k, t)*eval_mean(cumXY, k, t))/(eval_mean(cumSX, k, t) - (eval_mean(cumX, k, t))^2)
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t)*t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 + eval_mean(cumSY, k, t)) + costQ[shift(k-1)] + beta
  return(eval)
}

#############################################

evalReg_q <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1, t2) ### minimum of q_{t}^{k}, data y_{k} to y_{t}
{
  if(t == k){return(Inf)}
  eval <- (t - k + 1) * (eval_mean(cumSX, k, t) * t1^2 + 2 * eval_mean(cumX, k, t)*t1 * t2 + t2^2)
  eval <- eval -  (t - k + 1) * (2 * eval_mean(cumXY, k, t) * t1 + 2 * eval_mean(cumY, k, t) * t2 + eval_mean(cumSY, k, t)) + costQ[shift(k-1)] + beta
  return(eval)
}

##############################################
#####
##### Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + G #####
#####
##### s in [0,1]
#####
##############################################

ellipseBranch <- function(A, B, C, D, E, G, s, type)
{
  ###
  ### y in ]yMin, yMax[ to have two roots in x (for ellipses)
  ###
  u <- B^2 - A*C
  v <- B*D - A*E
  w <- v^2 - u*(D^2 - A*G)

  if((u > 0) | (w < 0)){return(NULL)}
  if((s < 0) | (s > 1)){return(NULL)}

  yMin <- (-v-sqrt(w))/u
  yMax <- (-v+sqrt(w))/u
  y <- yMin + s*(yMax - yMin)

  ###
  ### computing x value
  ###
  z1 <- B*y+D
  z2 <- z1^2 - A*(C*y^2 + 2*E*y + G)

  if(type == "left"){x <- (-z1 - sqrt(z2))/A}
  if(type == "right"){x <- (-z1 + sqrt(z2))/A}
  return(x)
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


