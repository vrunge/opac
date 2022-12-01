

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
  return(eval_mean(cumyS, j,k) - (eval_mean(cumy1,j,k)^2 + eval_mean(cumy2,j,k)^2))
}

#############################################

eval2D_q_min <- function(costQ, cumy1, cumy2, cumyS, k, t, beta) ###minimum of q_{t}^{k}, data y_k to y_t
{
  if(k == t){return(costQ[shift(k-1)] + beta)} ###costQ[shift(k)-1] = m_{k-1}
  return((t-k+1)*eval2D_var(cumy1, cumy2, cumyS, k,t) + costQ[shift(k-1)] + beta)
}

#############################################

eval2D_q <- function(costQ, cumy1, cumy2, cumyS, k, t, beta, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
{
  return((t - k + 1)*((t1 - eval_mean(cumy1, k, t))^2 +
                      (t2 - eval_mean(cumy2, k, t))^2) +
           eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta))
}

#############################################

eval2D_q_1_argmin <- function(R2, costQ, cumy1, cumy2, cumyS, j, k, t, beta) ###value of m_{t}^{jk}
{
  R <- sqrt(R2)
  t1kt <- eval_mean(cumy1, k, t)
  t2kt <- eval_mean(cumy2, k, t)
  t1jk <- eval_mean(cumy1, j, k-1)
  t2jk <- eval_mean(cumy2, j, k-1)

  D <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)
  return(eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) + (t - k + 1)*(D - R)^2)
}

#############################################

test2D_inclusion <- function(costQ, cumy1, cumy2, cumyS, j, k, t) #j < k to prune j
{
  Rsmall <- sqrt((costQ[shift(t-1)] - costQ[shift(j-1)])/(t-j) - eval2D_var(cumy1, cumy2, cumyS, j, t-1))
  Rbig <- sqrt((costQ[shift(t-1)] - costQ[shift(k-1)])/(t-k) - eval2D_var(cumy1, cumy2, cumyS, k, t-1))

  t1jt <- eval_mean(cumy1, j, t-1)
  t2jt <- eval_mean(cumy2, j, t-1)
  t1kt <- eval_mean(cumy1, k, t-1)
  t2kt <- eval_mean(cumy2, k, t-1)

  D <- sqrt((t1kt - t1jt)^2 + (t2kt - t2jt)^2)
  return((D + Rsmall) <= Rbig)
}


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#########
eval2D_q_0 <- function(costQ, cumy1, cumy2, cumyS, k, t, beta) ###minimum of q_{t}^{k}, data y_k to y_t
{
  if(k == t){return(list(p = c(eval_mean(cumy1, k, t),eval_mean(cumy2, k, t)),
                         m = costQ[shift(k-1)] + beta))} ###costQ[shift(k-1)] = m_{k-1}
  return(list(p = c(eval_mean(cumy1, k, t), eval_mean(cumy2, k, t)),
              m = (t - k + 1)*eval2D_var(cumy1, cumy2, cumyS, k, t) + costQ[shift(k-1)] + beta))
}


#########
eval2D_q_1 <- function(costQ, cumy1, cumy2, cumyS, j, k, t, beta) ###value of m_{t}^{jk}
{
  R2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval2D_var(cumy1, cumy2, cumyS, j, k-1)
  if(R2 < 0){return(list(m = Inf, p = c(Inf, Inf)))}
  R <- sqrt(R2)
  t1kt <- eval_mean(cumy1, k, t)
  t2kt <- eval_mean(cumy2, k, t)
  t1jk <- eval_mean(cumy1, j, k-1)
  t2jk <- eval_mean(cumy2, j, k-1)
  D <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)
  if(D == 0){return(list(p = c(t1kt + R, t2kt),
                         m = eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) + (t - k + 1)*R^2))}
  s <- R/D
  return(list(p = c(s*t1kt + (1-s)*t1jk, s*t2kt + (1-s)*t2jk),
              m =  eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) + (t - k + 1)*(D - R)^2))
}

#########
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
  if((D > Ri + Rj) | (D < Ri - Rj) | (D < Rj - Ri)){return(list(m = Inf, p = c(Inf, Inf)))}
  rho <- (1/2)*(Ri^2 - Rj^2)/D^2
  mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))
  t1A <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) + mu*(t2i - t2j)
  t2A <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) - mu*(t1i - t1j)
  mA <- eval2D_q(costQ, cumy1, cumy2, cumyS, k, t, beta, t1A, t2A)
  return(list(p = c(t1A, t2A), m = mA))
}

#########
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
  if((D > Ri + Rj) | (D < Ri - Rj) | (D < Rj - Ri)){return(list(m = Inf, p = c(Inf, Inf)))}
  rho <- (1/2)*(Ri^2 - Rj^2)/D^2
  mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))
  t1B <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) - mu*(t2i - t2j)
  t2B <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) + mu*(t1i - t1j)
  mB <- eval2D_q(costQ, cumy1, cumy2, cumyS, k, t, beta, t1B, t2B)
  return(list(p = c(t1B, t2B), m = mB))
}


#############################################
#############################################
############         Reg         ############
#############################################
#############################################






