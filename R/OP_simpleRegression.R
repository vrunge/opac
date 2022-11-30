
############################################
#############     OP_Reg    ################
############################################

#' OP_Reg
#' @description Optimal Partitioning algorithm for changes in simple regression (no pruning)
#' @param data a dataframe with two components: x and y, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment)
#' @examples
#' OP_Reg(dataGenerator_Reg(chpts = c(50,100,150), A = c(-1,1,-1), B = c(-1,1,-1)))
OP_Reg <- function(data, beta = 4 * log(nrow(data)))
{
  #########
  ###
  ### DATA preprocessing
  ###
  x <- data$x
  y <- data$y
  n <- length(y)
  cumX <- cumsum(c(0, x))
  cumY <- cumsum(c(0, y))
  cumXY <- cumsum(c(0, x*y))
  cumSX <- cumsum(c(0, x^2))
  cumSY <- cumsum(c(0, y^2))

  #########
  ###
  ### INNER FUNCTIONS (EVAL)
  ###
  shift <- function(k){return(k+1)}  ###costQ(shift(k)) = m_k
  eval_meanX <- function(j, k){return((cumX[k+1]-cumX[j])/(k-j+1))}
  eval_meanY <- function(j, k){return((cumY[k+1]-cumY[j])/(k-j+1))}
  eval_meanX2 <- function(j, k){return((cumSX[k+1]-cumSX[j])/(k-j+1))}
  eval_meanY2 <- function(j, k){return((cumSY[k+1]-cumSY[j])/(k-j+1))}
  eval_meanXY <- function(j, k){return((cumXY[k+1]-cumXY[j])/(k-j+1))}

  #########
  eval_q_min <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(t == k){return(Inf)}
    t1 <- (eval_meanXY(k,t) - eval_meanX(k,t)*eval_meanY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    t2 <- (eval_meanX2(k,t)*eval_meanY(k,t) - eval_meanX(k,t)*eval_meanXY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    eval <-  (t-k+1)*eval_meanX2(k,t)*t1^2 + 2*(t-k+1)*eval_meanX(k,t)*t1*t2 + (t-k+1)*t2^2
    eval <- eval - 2*(t-k+1)*eval_meanY(k,t)*t2 - 2*(t-k+1)*eval_meanXY(k,t)*t1 + (t-k+1)*eval_meanY2(k,t) + costQ[shift(k)-1] + beta
    return(eval)
  }

  #########
  ###
  ### INITIALIZATION
  ###
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  index <- 0 #best index (to be found) for last change-point

  #########
  ###
  ### UPDATE rule Dynamic Programming
  ###
  index <- 1
  for(t in 1:n)
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
  return(list(changepoints = changepoints[-1], nb = NULL))
}





#################################################
#############     OP_Reg_PELT    ################
#################################################

#' OP_Reg_PELT
#' @description Optimal Partitioning algorithm for changes in simple regression (with PELT pruning)
#' @param data a dataframe with two components: x and y, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_Reg_PELT(dataGenerator_Reg(chpts = c(50,100,150), A = c(-1,1,-1), B = c(-1,1,-1)))
OP_Reg_PELT <- function(data, beta = 4 * log(nrow(data)))
{
  #########
  ###
  ### DATA preprocessing
  ###
  x <- data$x
  y <- data$y
  n <- length(y)
  cumX <- cumsum(c(0, x))
  cumY <- cumsum(c(0, y))
  cumXY <- cumsum(c(0, x*y))
  cumSX <- cumsum(c(0, x^2))
  cumSY <- cumsum(c(0, y^2))

  #########
  ###
  ### INNER FUNCTIONS (EVAL)
  ###
  shift <- function(k){return(k+1)}  ###costQ(shift(k)) = m_k
  eval_meanX <- function(j, k){return((cumX[k+1]-cumX[j])/(k-j+1))}
  eval_meanY <- function(j, k){return((cumY[k+1]-cumY[j])/(k-j+1))}
  eval_meanX2 <- function(j, k){return((cumSX[k+1]-cumSX[j])/(k-j+1))}
  eval_meanY2 <- function(j, k){return((cumSY[k+1]-cumSY[j])/(k-j+1))}
  eval_meanXY <- function(j, k){return((cumXY[k+1]-cumXY[j])/(k-j+1))}

  #########
  eval_q_min <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(t == k){return(Inf)}
    t1 <- (eval_meanXY(k,t) - eval_meanX(k,t)*eval_meanY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    t2 <- (eval_meanX2(k,t)*eval_meanY(k,t) - eval_meanX(k,t)*eval_meanXY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    eval <-  (t-k+1)*eval_meanX2(k,t)*t1^2 + 2*(t-k+1)*eval_meanX(k,t)*t1*t2 + (t-k+1)*t2^2
    eval <- eval - 2*(t-k+1)*eval_meanY(k,t)*t2 - 2*(t-k+1)*eval_meanXY(k,t)*t1 + (t-k+1)*eval_meanY2(k,t) + costQ[k] + beta
    return(eval)
  }

  eval_segment_kt <- function(k, t)
  {
    if(t == k){return(0)}
    t1 <- (eval_meanXY(k,t) - eval_meanX(k,t)*eval_meanY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    t2 <- (eval_meanX2(k,t)*eval_meanY(k,t) - eval_meanX(k,t)*eval_meanXY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    eval <-  (t-k+1)*eval_meanX2(k,t)*t1^2 + 2*(t-k+1)*eval_meanX(k,t)*t1*t2 + (t-k+1)*t2^2
    eval <- eval - 2*(t-k+1)*eval_meanY(k,t)*t2 - 2*(t-k+1)*eval_meanXY(k,t)*t1 + (t-k+1)*eval_meanY2(k,t)
    return(eval)
  }

  #########
  ###
  ### INITIALIZATION
  ###
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- 1

  #########
  ###
  ### UPDATE rule Dynamic Programming
  ###
  for(t in 1:n)
  {
    min_temp <- Inf
    index <- 1
    for(k in indexSet)
    {
      eval <- eval_segment_kt(k,t) + costQ[k] + beta
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1

    ### PRUNING STEP
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      if(costQ[shift(t)]> costQ[k] + eval_segment_kt(k,t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    indexSet <- c(nonpruned, t+1) #add new test point
    nb[t] <- length(indexSet) ### count number of elements
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




###############################################
#############     OP_Reg_1C    ################
###############################################

#' OP_Reg_1C
#' @description Optimal Partitioning algorithm for changes in simple regression (with OP1C pruning)
#' @param data a dataframe with two components: x and y, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_Reg_1C(dataGenerator_Reg(chpts = c(50,100,150), A = c(-1,1,-1), B = c(-1,1,-1)))
OP_Reg_1C <- function(data, beta = 4 * log(nrow(data)))
{
  #########
  ###
  ### DATA preprocessing
  ###
  x <- data$x
  y <- data$y
  n <- length(y)
  cumX <- cumsum(c(0, x))
  cumY <- cumsum(c(0, y))
  cumXY <- cumsum(c(0, x*y))
  cumSX <- cumsum(c(0, x^2))
  cumSY <- cumsum(c(0, y^2))

  #########
  ###
  ### INNER FUNCTIONS (EVAL)
  ###
  shift <- function(k){return(k+1)}  ###costQ(shift(k)) = m_k
  eval_meanX <- function(j, k){return((cumX[k+1]-cumX[j])/(k-j+1))}
  eval_meanY <- function(j, k){return((cumY[k+1]-cumY[j])/(k-j+1))}
  eval_meanX2 <- function(j, k){return((cumSX[k+1]-cumSX[j])/(k-j+1))}
  eval_meanY2 <- function(j, k){return((cumSY[k+1]-cumSY[j])/(k-j+1))}
  eval_meanXY <- function(j, k){return((cumXY[k+1]-cumXY[j])/(k-j+1))}

  #########
  eval_q_min <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(t == k){return(Inf)}
    t1 <- (eval_meanXY(k,t) - eval_meanX(k,t)*eval_meanY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    t2 <- (eval_meanX2(k,t)*eval_meanY(k,t) - eval_meanX(k,t)*eval_meanXY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    eval <-  (t-k+1)*eval_meanX2(k,t)*t1^2 + 2*(t-k+1)*eval_meanX(k,t)*t1*t2 + (t-k+1)*t2^2
    eval <- eval - 2*(t-k+1)*eval_meanY(k,t)*t2 - 2*(t-k+1)*eval_meanXY(k,t)*t1 + (t-k+1)*eval_meanY2(k,t) + costQ[k] + beta
    return(eval)
  }

  eval_segment_kt <- function(k, t)
  {
    if(t == k){return(0)}
    t1 <- (eval_meanXY(k,t) - eval_meanX(k,t)*eval_meanY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    t2 <- (eval_meanX2(k,t)*eval_meanY(k,t) - eval_meanX(k,t)*eval_meanXY(k,t))/(eval_meanX2(k,t) - (eval_meanX(k,t))^2)
    eval <-  (t-k+1)*eval_meanX2(k,t)*t1^2 + 2*(t-k+1)*eval_meanX(k,t)*t1*t2 + (t-k+1)*t2^2
    eval <- eval - 2*(t-k+1)*eval_meanY(k,t)*t2 - 2*(t-k+1)*eval_meanXY(k,t)*t1 + (t-k+1)*eval_meanY2(k,t)
    return(eval)
  }

  #########
  ###
  ### INITIALIZATION
  ###
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- 1

  #########
  ###
  ### UPDATE rule Dynamic Programming
  ###
  for(t in 1:n)
  {
    min_temp <- Inf
    index <- 1
    for(k in indexSet)
    {
      eval <- eval_segment_kt(k,t) + costQ[k] + beta
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1

    ### PRUNING STEP
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      if(costQ[shift(t)]> costQ[k] + eval_segment_kt(k,t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    indexSet <- c(nonpruned, t+1) #add new test point
    nb[t] <- length(indexSet) ### count number of elements
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


