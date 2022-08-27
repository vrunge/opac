

############################################
#############     OP_2D    #################
############################################

#' OP_2D
#' @description Optimal Partitioning algorithm for bivariate independent time series (no pruning)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment)
#' @examples
#' OP_2D(dataGenerator2D(chpts = c(30,100,120), means1 = c(0,5,0), means2 = c(7,1,-4)))
OP_2D <- function(data, beta = 4 * log(nrow(data)))
{
  #########
  ###
  ### DATA preprocessing
  ###
  y1 <- data$y1
  y2 <- data$y2
  n <- length(y1)
  cumy1 <- cumsum(c(0, y1))
  cumy2 <- cumsum(c(0, y2))
  cumyS <- cumsum(c(0, y1^2 + y2^2))

  #########
  ###
  ### INNER FUNCTIONS (EVAL)
  ###
  shift <- function(k){return(k+1)}  ###costQ(shift(k)) = m_k
  eval_meany1 <- function(j, k){return((cumy1[k+1]-cumy1[j])/(k-j+1))}
  eval_meany2 <- function(j, k){return((cumy2[k+1]-cumy2[j])/(k-j+1))}
  eval_meanS <- function(j, k){return((cumyS[k+1]-cumyS[j])/(k-j+1))}

  #########
  eval_var <- function(j, k) ###Variance for datapoint y_j to y_{k}
  {
    if(j == k){return(0)}
    return(eval_meanS(j,k) - (eval_meany1(j,k)^2 + eval_meany2(j,k)^2))
  }
  #########
  eval_q_min <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(k == t){return(costQ[shift(k)-1] + beta)} ###costQ[shift(k)-1] = m_{k-1}
    return((shift(t)-k)*eval_var(k,t) + costQ[shift(k)-1] + beta)
  }

  #########
  ###
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  index <- 0 #best index (to be found) for last change-point

  #########
  ###
  ### UPDATE rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    min_temp <- Inf
    for(k in 1:t)
    {
      eval <- eval_q_min(k,t)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1
  }

  #########
  ###
  ### backtracking step
  ###
  cp <- cp[-1] # remove first value
  changepoints <- n # vector of change-point to build
  current <- n

  while(current > 1)
  {
    pointval <- cp[current] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1]))
}


#############################################
############     OP_2D_PELT    ##############
#############################################

#' OP_2D_PELT
#' @description Optimal Partitioning algorithm for bivariate independent time series (with PELT pruning)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_PELT(dataGenerator2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_PELT <- function(data, beta = 4 * log(nrow(data)))
{
  #########
  ###
  ### DATA preprocessing
  ###
  y1 <- data$y1
  y2 <- data$y2
  n <- length(y1)
  cumy1 <- cumsum(c(0, y1))
  cumy2 <- cumsum(c(0, y2))
  cumyS <- cumsum(c(0, y1^2 + y2^2))

  #########
  ###
  ### INNER FUNCTIONS (EVAL)
  ###
  shift <- function(k){return(k+1)}  ###costQ(shift(k)) = m_k
  eval_meany1 <- function(j, k){return((cumy1[k+1]-cumy1[j])/(k-j+1))}
  eval_meany2 <- function(j, k){return((cumy2[k+1]-cumy2[j])/(k-j+1))}
  eval_meanS <- function(j, k){return((cumyS[k+1]-cumyS[j])/(k-j+1))}

  #########
  eval_var <- function(j, k) ###Variance for datapoint y_j to y_{k}
  {
    if(j == k){return(0)}
    return(eval_meanS(j,k) - (eval_meany1(j,k)^2 + eval_meany2(j,k)^2))
  }
  #########
  eval_q_min <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(k == t){return(costQ[shift(k)-1] + beta)} ###costQ[shift(k)-1] = m_{k-1}
    return((t-k+1)*eval_var(k,t) + costQ[shift(k)-1] + beta)
  }

  #########
  ###
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta

  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- 1

  #########
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- eval_q_min(k,t)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1

    ### PRUNING STEP
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      if(costQ[shift(t)] > costQ[shift(k-1)] + (t-k+1)*eval_var(k,t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    nb[t] <- length(indexSet) ### count number of elements
    indexSet <- c(nonpruned, shift(t)) #add new test point
  }

  #########
  ###
  ### backtracking step
  ###
  cp <- cp[-1] # remove first value
  changepoints <- n # vector of change-point to build
  current <- n

  while(current > 1)
  {
    pointval <- cp[current] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1], nb = nb[1:(n-1)]))
}
