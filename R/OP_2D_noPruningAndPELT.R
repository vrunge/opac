

############################################
#############     OP_2D    #################
############################################

#' OP_2D
#' @description Optimal Partitioning algorithm for bivariate independent time series (no pruning)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment)
#' @examples
#' OP_2D(dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,5,0), means2 = c(7,1,-4)))
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
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  index <- 0 #best index (to be found) for last change-point

  #########
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    min_temp <- Inf

    for(k in 1:t)
    {
      eval <- eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta)
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

#############################################
############     OP_2D_PELT    ##############
#############################################

#' OP_2D_PELT
#' @description Optimal Partitioning algorithm for bivariate independent time series (with PELT pruning)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_PELT(dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
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
  ### INITIALIZATION
  ###
  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta

  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- NULL

  #########
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    indexSet <- c(indexSet, t) #add new test point

    min_temp <- Inf

    for(k in indexSet)
    {
      eval <- eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t)] <- min_temp
    cp[shift(t)] <- index - 1


    ### PRUNING STEP
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inequality-based)
      if(costQ[shift(t)] >= costQ[shift(k-1)] + (t-k+1)*eval2D_var(cumy1, cumy2, cumyS, k, t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    indexSet <- nonpruned
    nb[t] <- length(indexSet) ### count number of elements
    #print("zzzzzzzzzzzz")
    #print(t)
    #print(indexSet)
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
