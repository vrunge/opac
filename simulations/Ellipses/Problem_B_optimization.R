
#y1 <- data$y1
#y2 <- data$y2
#n <- length(y1)
#cumy1 <- cumsum(c(0, y1))
#cumy2 <- cumsum(c(0, y2))
#cumyS <- cumsum(c(0, y1^2 + y2^2))

#costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
#costQ[1] <- -beta #costQ[2] = 0

shift <- function(k){return(k + 1)}  ###costQ(shift(k)) = m_k
evalB_meany1 <- function(cumy1, j, k){return((cumy1[k+1] - cumy1[j])/(k-j+1))}
evalB_meany2 <- function(cumy2, j, k){return((cumy2[k+1] - cumy2[j])/(k-j+1))}
evalB_meanSq <- function(cumyS, j, k){return((cumyS[k+1] - cumyS[j])/(k-j+1))}

evalB_var <- function(cumy1, cumy2, cumyS, j, k) ###Variance for datapoint y_j to y_k
{
  if(j == k){return(0)}
  return(evalB_meanSq(cumyS, j, k) - (evalB_meany1(cumy1, j, k)^2 + evalB_meany2(cumy2, j, k)^2))
}

### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###

evalB_q <- function(cumy1, cumy2, cumyS, costQ, k, t, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
{
  return((t - k + 1)*((t1 - evalB_meany1(cumy1, k, t))^2 +  (t2 - evalB_meany2(cumy2, k, t))^2)
         + evalB_q(cumy1, cumy2, cumyS, costQ, k, t)$m)
}

### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ###

evalB_q_0 <- function(cumy1, cumy2, cumyS, costQ, k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
{
  return(list(p = c(0,0),
              m =  0))
}


#########
evalB_q_1 <- function(cumy1, cumy2, cumyS, costQ, j, k, t) ###value of m_{t}^{jk}
{

  return(list(p = c(0,0),
              m =  0))
}

#########
evalB_q_2 <- function(cumy1, cumy2, cumyS, costQ, i, j, k, t) ###value of m_{t}^{ijk}
{
  return(list(p = c(0,0),
              m =  0))
}
