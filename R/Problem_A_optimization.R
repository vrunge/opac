
#y1 <- data$y1
#y2 <- data$y2
#n <- length(y1)
#cumy1 <- cumsum(c(0, y1))
#cumy2 <- cumsum(c(0, y2))
#cumyS <- cumsum(c(0, y1^2 + y2^2))

#costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
#costQ[1] <- -beta #costQ[2] = 0

shift <- function(k){return(k + 1)}  ###costQ(shift(k)) = m_k
evalA_meany1 <- function(cumy1, j, k){return((cumy1[k+1] - cumy1[j])/(k-j+1))}
evalA_meany2 <- function(cumy2, j, k){return((cumy2[k+1] - cumy2[j])/(k-j+1))}
evalA_meanSq <- function(cumyS, j, k){return((cumyS[k+1] - cumyS[j])/(k-j+1))}

evalA_var <- function(cumy1, cumy2, cumyS, j, k) ###Variance for datapoint y_j to y_k
{
  if(j == k){return(0)}
  return(evalA_meanSq(cumyS, j, k) - (evalA_meany1(cumy1, j, k)^2 + evalA_meany2(cumy2, j, k)^2))
}

### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###

evalA_q <- function(cumy1, cumy2, cumyS, costQ, k, t, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
{
  return((t - k + 1)*((t1 - evalA_meany1(cumy1, k, t))^2 +  (t2 - evalA_meany2(cumy2, k, t))^2)
         + evalA_q(cumy1, cumy2, cumyS, costQ, k, t)$m)
}

### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###

evalA_q_0 <- function(cumy1, cumy2, cumyS, costQ, k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
{
  if(k == t){return(list(p = c(evalA_meany1(cumy1, t, t), evalA_meany2(cumy2, t, t)),
                         m = costQ[shift(t-1)] + beta))} ###costQ[shift(k)-1] = m_{k-1}

  return(list(p = c(evalA_meany1(cumy1,k,t), evalA_meany2(cumy2, k, t)),
              m = (t - k + 1)*evalA_var(cumy1, cumy2, cumyS, k, t) + costQ[shift(k-1)] + beta))
}

#########
evalA_q_1 <- function(cumy1, cumy2, cumyS, costQ, j, k, t) ###value of m_{t}^{jk}
{
  R <- sqrt((costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - evalA_var(cumy1, cumy2, cumyS, j, k-1))
  t1kt <- evalA_meany1(cumy1, k, t)
  t2kt <- evalA_meany2(cumy2, k, t)
  t1jk <- evalA_meany1(cumy1, j, k-1)
  t2jk <- evalA_meany2(cumy2, j, k-1)
  D <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)

  if(D == 0){return(list(p = c(t1kt + R, t2kt),
                         m = evalA_q_0(cumy1, cumy2, cumyS, costQ, k, t)$m + (t - k + 1)*R^2))}
  s <- R/D
  return(list(p = c(s*t1kt + (1-s)*t1jk, s*t2kt + (1-s)*t2jk),
              m =  evalA_q_0(cumy1, cumy2, cumyS, costQ, k, t)$m + (t - k + 1)*(D - R)^2))
}

#########
evalA_q_2 <- function(cumy1, cumy2, cumyS, costQ, i, j, k, t) ###value of m_{t}^{ijk}
{
  Ri2 <- (costQ[shift(k-1)] - costQ[shift(i-1)])/(k-i) - evalA_var(cumy1, cumy2, cumyS, i, k-1)
  Rj2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - evalA_var(cumy1, cumy2, cumyS, j, k-1)
  if((Ri2 < 0) || (Rj2 < 0)){return(list(m = Inf, p = c(0,0)))}
  Ri <- sqrt(Ri2)
  Rj <- sqrt(Rj2)
  t1i <- evalA_meany1(cumy1, i, k-1)
  t2i <- evalA_meany2(cumy2, i, k-1)
  t1j <- evalA_meany1(cumy1, j, k-1)
  t2j <- evalA_meany2(cumy2, j, k-1)
  D <- sqrt((t1i - t1j)^2 + (t2i - t2j)^2)

  if((D > Ri + Rj) | (D < Ri - Rj) | (D < Rj - Ri))
    {return(list(p1 = c(0,0), m1 = Inf, p2 = c(0,0), m2 = Inf))}

  rho <- (1/2)*(Ri^2 - Rj^2)/D^2
  mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))
  t1A <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) + mu*(t2i - t2j)
  t2A <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) - mu*(t1i - t1j)

  t1B <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) - mu*(t2i - t2j)
  t2B <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) + mu*(t1i - t1j)

  mA <- evalA_q(cumy1, cumy2, cumyS, costQ, k, t, t1A, t2A)
  mB <- evalA_q(cumy1, cumy2, cumyS, costQ, k, t, t1B, t2B)

  return(list(p1 = c(t1A, t2A), m1 = mA, p2 = c(t1B, t2B), m2 = mB))
}
