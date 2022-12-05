
############################################
#############     OP_Reg    ################
############################################

#' OP_Reg
#'
#' @description Optimal Partitioning algorithm for changes in simple regression (no pruning)
#' @param data a data-frame with two components: x and y, time series of same length
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
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  index <- 1 #best index (to be found) for last change-point

  #########
  ###
  ### UPDATE rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    min_temp <- Inf

    for(k in 1:t)
    {
      eval <- evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1
  }

  ######### DDDDDDDDDDDDDDDDDDDDDD #########
  ###
  ### backtracking step
  ###
  changepoints <- n # vector of change-point to build
  current <- n

  while(changepoints[1] > 0)
  {
    pointval <- cp[shift(current)] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1]))
}


#################################################
#############     OP_Reg_PELT    ################
#################################################

#' OP_Reg_PELT
#' @description Optimal Partitioning algorithm for changes in simple regression (with PELT pruning)
#' @param data a data-frame with two components: x and y, time series of same length
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
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
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
      eval <- evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1


    ### PRUNING STEP
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      ### strict ">" is better (for pruning)
      if(costQ[shift(t)] >= costQ[shift(k-1)] + evalReg_segment_kt(cumX, cumY, cumXY, cumSX, cumSY, k, t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    indexSet <- c(nonpruned, t+1) #add new point
    nb[t] <- length(indexSet) ### count number of elements
  }

  ######### DDDDDDDDDDDDDDDDDDDDDD #########
  ###
  ### backtracking step
  ###
  changepoints <- n # vector of change-point to build
  current <- n

  while(changepoints[1] > 0)
  {
    pointval <- cp[shift(current)] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1],  nb = nb))
}

