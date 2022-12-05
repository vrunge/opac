

###############################################################
####################### INNER FUNCTIONS #######################
###############################################################

shift <- function(k){return(k+1)}
eval_mean <- function(v, j, k){return((v[k+1]-v[j])/(k-j+1))} ##mean from v_j to v_k included


#############################################
#############################################
############         2D          ############
#############################################
#############################################


eval2D_var <- function(cumy1, cumy2, cumyS, j, k) ###Variance for data-point y_j to y_k (over the 2 dimensions)
{
  if(j == k){return(0)}
  return(eval_mean(cumyS, j, k) - (eval_mean(cumy1, j, k)^2 + eval_mean(cumy2, j, k)^2))
}

#############################################

eval2D_q_min <- function(costQ, cumy1, cumy2, cumyS, k, t, beta) ###minimum of q_{t}^{k}, data y_k to y_t
{
  if(k == t){return(costQ[shift(k-1)] + beta)} ###costQ[shift(k)-1] = m_{k-1}
  return((t-k+1)*eval2D_var(cumy1, cumy2, cumyS, k, t) + costQ[shift(k-1)] + beta)
}

#############################################

eval2D_q <- function(costQ, cumy1, cumy2, cumyS, k, t, beta, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
{
  return((t - k + 1)*((t1 - eval_mean(cumy1, k, t))^2 +
                      (t2 - eval_mean(cumy2, k, t))^2) +
           eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta))
}

#############################################

eval2D_q_1_min <- function(R2, costQ, cumy1, cumy2, cumyS, j, k, t, beta) ###value of m_{t}^{jk}
{
  R <- sqrt(R2)
  t1jt <- eval_mean(cumy1, j, t)
  t2jt <- eval_mean(cumy2, j, t)
  t1kt <- eval_mean(cumy1, k, t)
  t2kt <- eval_mean(cumy2, k, t)
  t1jk <- eval_mean(cumy1, j, k-1)
  t2jk <- eval_mean(cumy2, j, k-1)

  Dj <- sqrt((t1jt - t1jk)^2 + (t2jt - t2jk)^2)
  Dk <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)

  mj <- eval2D_q_min(costQ, cumy1, cumy2, cumyS, j, t, beta) + (t - j + 1)*(Dj - R)^2
  mk <- eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) + (t - k + 1)*(Dk - R)^2

  res <- min(mj, mk)
  res <- trunc(res*1e13)*1e-13
  return(res)
}


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#############################################

eval2D_q_0 <- function(costQ, cumy1, cumy2, cumyS, k, t, beta) ###minimum of q_{t}^{k}, data y_k to y_t
{
  if(k == t){return(list(p = c(eval_mean(cumy1, k, t),eval_mean(cumy2, k, t)),
                         m = costQ[shift(k-1)] + beta))} ###costQ[shift(k-1)] = m_{k-1}
  return(list(p = c(eval_mean(cumy1, k, t), eval_mean(cumy2, k, t)),
              m = (t - k + 1)*eval2D_var(cumy1, cumy2, cumyS, k, t) + costQ[shift(k-1)] + beta))
}

#############################################

eval2D_q_1 <- function(costQ, cumy1, cumy2, cumyS, j, k, t, beta) ###value of m_{t}^{jk}
{
  R2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval2D_var(cumy1, cumy2, cumyS, j, k-1)
  if(R2 < 0){return(list(m = Inf, p = c(Inf, Inf)))}
  R <- sqrt(R2)
  t1jt <- eval_mean(cumy1, j, t)
  t2jt <- eval_mean(cumy2, j, t)
  t1kt <- eval_mean(cumy1, k, t)
  t2kt <- eval_mean(cumy2, k, t)
  t1jk <- eval_mean(cumy1, j, k-1)
  t2jk <- eval_mean(cumy2, j, k-1)
  Dj <- sqrt((t1jt - t1jk)^2 + (t2jt - t2jk)^2)
  Dk <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)

  if(Dj == 0){return(list(p = c(t1jt + R, t2jt),
                          m = eval2D_q_min(costQ, cumy1, cumy2, cumyS, j, t, beta) + (t - j + 1)*R^2))}
  if(Dk == 0){return(list(p = c(t1kt + R, t2kt),
                          m = eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) + (t - k + 1)*R^2))}

  mj <- eval2D_q_min(costQ, cumy1, cumy2, cumyS, j, t, beta) + (t - j + 1)*(Dj - R)^2
  mk <- eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) + (t - k + 1)*(Dk - R)^2
  choix <- which.min(c(mj, mk))
  if(choix == 1)
  {
    s <- R/Dj
    return(list(p = c(s*t1jt + (1-s)*t1jk, s*t2jt + (1-s)*t2jk), m =  trunc(mj*1e13)*1e-13))
  }
  if(choix == 2)
  {
    s <- R/Dk
    return(list(p = c(s*t1kt + (1-s)*t1jk, s*t2kt + (1-s)*t2jk), m =  trunc(mk*1e13)*1e-13))
  }
}

#############################################

eval2D_q_21 <- function(costQ, cumy1, cumy2, cumyS, i, j, k, t, beta) ###value of m_{t}^{ijk}
{
  Ri2 <- (costQ[shift(k-1)] - costQ[shift(i-1)])/(k-i) - eval2D_var(cumy1, cumy2, cumyS, i, k-1)
  Rj2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval2D_var(cumy1, cumy2, cumyS, j, k-1)
  if((Ri2 < 0) || (Rj2 < 0)){return(list(m = Inf, p = c(Inf, Inf)))}
  Ri <- sqrt(Ri2)
  Rj <- sqrt(Rj2)
  t1i <- eval_mean(cumy1, i, k-1)
  t2i <- eval_mean(cumy2, i, k-1)
  t1j <- eval_mean(cumy1, j, k-1)
  t2j <- eval_mean(cumy2, j, k-1)
  D <- sqrt((t1i - t1j)^2 + (t2i - t2j)^2)

  if((D == 0) | (D > Ri + Rj) | (D < abs(Ri - Rj))){return(list(m = Inf, p = c(Inf, Inf)))}

  rho <- (1/2)*(Rj^2 - Ri^2)/D^2
  mu <- (1/(2*D^2))*sqrt((Rj^2 - (D + Ri)^2)*((Ri - D)^2 - Rj^2))
  #mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))

  t1A <- (1/2)*(t1i + t1j) + rho*(t1i - t1j) + mu*(t2j - t2i)
  t2A <- (1/2)*(t2i + t2j) + rho*(t2i - t2j) - mu*(t1j - t1i)

  mAi <- eval2D_q(costQ, cumy1, cumy2, cumyS, i, t, beta, t1A, t2A)
  mAj <- eval2D_q(costQ, cumy1, cumy2, cumyS, j, t, beta, t1A, t2A)
  mAk <- eval2D_q(costQ, cumy1, cumy2, cumyS, k, t, beta, t1A, t2A)

  res <- min(mAi, mAj, mAk)
  res <- trunc(res*1e13)*1e-13

  return(list(p = c(t1A, t2A), m = res))
}

#############################################

eval2D_q_22 <- function(costQ, cumy1, cumy2, cumyS, i, j, k, t, beta) ###value of m_{t}^{ijk}
{
  Ri2 <- (costQ[shift(k-1)] - costQ[shift(i-1)])/(k-i) - eval2D_var(cumy1, cumy2, cumyS, i, k-1)
  Rj2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval2D_var(cumy1, cumy2, cumyS, j, k-1)
  if((Ri2 < 0) || (Rj2 < 0)){return(list(m = Inf, p = c(Inf, Inf)))}
  Ri <- sqrt(Ri2)
  Rj <- sqrt(Rj2)
  t1i <- eval_mean(cumy1, i, k-1)
  t2i <- eval_mean(cumy2, i, k-1)
  t1j <- eval_mean(cumy1, j, k-1)
  t2j <- eval_mean(cumy2, j, k-1)
  D <- sqrt((t1i - t1j)^2 + (t2i - t2j)^2)

  if((D == 0) | (D > Ri + Rj) | (D < abs(Ri - Rj))){return(list(m = Inf, p = c(Inf, Inf)))}

  rho <- (1/2)*(Rj^2 - Ri^2)/D^2
  mu <- (1/(2*D^2))*sqrt((Rj^2 - (D + Ri)^2)*((Ri - D)^2 - Rj^2))
  #mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))
  t1B <- (1/2)*(t1i + t1j) + rho*(t1i - t1j) - mu*(t2j - t2i)
  t2B <- (1/2)*(t2i + t2j) + rho*(t2i - t2j) + mu*(t1j - t1i)

  mBi <- eval2D_q(costQ, cumy1, cumy2, cumyS, i, t, beta, t1B, t2B)
  mBj <- eval2D_q(costQ, cumy1, cumy2, cumyS, j, t, beta, t1B, t2B)
  mBk <- eval2D_q(costQ, cumy1, cumy2, cumyS, k, t, beta, t1B, t2B)

  res <- min(mBi, mBj, mBk)
  res <- trunc(res*1e13)*1e-13

  return(list(p = c(t1B, t2B), m = res))
}

#############################################

eval2D_circleIntersection_ijk <- function(costQ, cumy1, cumy2, cumyS, i, j, k)
{
  Ri2 <- (costQ[shift(k-1)] - costQ[shift(i-1)])/(k-i) - eval2D_var(cumy1, cumy2, cumyS, i, k-1)
  Rj2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval2D_var(cumy1, cumy2, cumyS, j, k-1)
  if((Ri2 < 0) || (Rj2 < 0)){return(FALSE)}
  Ri <- sqrt(Ri2)
  Rj <- sqrt(Rj2)
  t1i <- eval_mean(cumy1, i, k-1)
  t2i <- eval_mean(cumy2, i, k-1)
  t1j <- eval_mean(cumy1, j, k-1)
  t2j <- eval_mean(cumy2, j, k-1)
  D <- sqrt((t1i - t1j)^2 + (t2i - t2j)^2)
  if((D > Ri + Rj) | (D < abs(Ri - Rj))){return(FALSE)}
  return(TRUE)
}


