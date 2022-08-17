


############################################
#############     OP_2D    #################
############################################

#' OP_2D
#' @description Optimal Partitioning algorithm for bivariate independent time series (no pruning)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty
#' @return a list with the change-point elements (each last index of each segment)
#' @examples
#' OP_2D(dataGenerator2D(chpts = c(30,100,120), means1 = c(0,5,0), means2 = c(7,1,-4)))
OP_2D <- function(data, beta = 4*log(length(data$y1)))
{
  ###
  ### INITIALIZATION
  ###
  y1 <- data$y1
  y2 <- data$y2
  n <- length(y1)
  cumy1 <- cumsum(c(0, y1))
  cumy2 <- cumsum(c(0, y2))
  cumyS <- cumsum(c(0, y1^2 + y2^2))

  cp <- rep(0, n + 1) #cp vector cp[i] = index of the last change-point for data y(1) to y(i-1)
  costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta
  index <- 0 #best index (to be found) for last change-point

  ###
  ### UPDATE rule Dynamic programming
  ###
  for(t in 1:n)
  {
    min_temp <- Inf
    for(i in 1:t)
    {
      val <- cumyS[t+1]-cumyS[i] - ((cumy1[t+1]-cumy1[i])^2+(cumy2[t+1]-cumy2[i])^2)/(t-i+1) + costQ[i] + beta
      if(val < min_temp)
      {
        min_temp <- val
        index <- i
      }
      costQ[t+1] <- min_temp
      cp[t+1] <- index - 1
    }
  }

  cp <- cp[-1] # remove first value
  ###
  ### backtracking step
  ###
  changepoints <- n # vector of change-point to build
  current <- n

  while(current > 1)
  {
    pointval <- cp[current] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1]))
  #return(list(cp = cp, costQ = costQ, changepoints = changepoints[-1]))
}


#############################################
############     OP_2D_PELT    ##############
#############################################

#' OP_2D_PELT
#' @description Optimal Partitioning algorithm for bivariate independent time series (with PELT pruning)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_PELT(dataGenerator2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_PELT <- function(data, beta = 4*log(length(data$y1)))
{
  ###
  ### INITIALIZATION
  ###
  y1 <- data$y1
  y2 <- data$y2
  n <- length(y1)
  cumy1 <- cumsum(c(0, y1))
  cumy2 <- cumsum(c(0, y2))
  cumyS <- cumsum(c(0, y1^2 + y2^2))

  cp <- rep(0, n + 1) #cp vector cp[i] = index of the last change-point for data y(1) to y(i-1)
  costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta
  index <- 0 #best index (to be found) for last change-point

  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- 1
  ###
  ### update rule Dynamic programming
  ###
  for(t in 1:n)
  {
    min_temp <- Inf
    for(i in indexSet)
    {
      val <- cumyS[t+1]-cumyS[i] - ((cumy1[t+1]-cumy1[i])^2+(cumy2[t+1]-cumy2[i])^2)/(t-i+1) + costQ[i] + beta
      if(val < min_temp)
      {
        min_temp <- val
        index <- i
      }
      costQ[t+1] <- min_temp
      cp[t+1] <- index - 1
    }

    ### PRUNING STEP
    nonpruned <- NULL
    for(i in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      if(costQ[t+1] > costQ[i] + cumyS[t+1]-cumyS[i] - ((cumy1[t+1]-cumy1[i])^2+(cumy2[t+1]-cumy2[i])^2)/(t-i+1))
      {
        nonpruned <- c(nonpruned, i)
      }
    }
    indexSet <- c(nonpruned, t+1) #add new test point
    nb[t] <- length(indexSet) ### count number of elements
  }

  cp <- cp[-1] # remove first value

  ###
  ### backtracking step
  ###
  changepoints <- n # vector of change-point to build
  current <- n

  while(current > 1)
  {
    pointval <- cp[current] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  #return(list(cp = cp, costQ = costQ, changepoints = changepoints[-1]))
  return(list(changepoints = changepoints[-1], nb = nb))
}



###############################################
#############     OP_2D_1C    #################
###############################################

#' OP_2D_1C
#' @description Optimal Partitioning algorithm for bivariate independent time series (with OP1C algorithm)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_1C(dataGenerator2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_1C <- function(data, beta = 4*log(length(data$y1)))
{
  ###
  ### INITIALIZATION
  ###
  y1 <- data$y1
  y2 <- data$y2
  n <- length(y1)
  cumy1 <- cumsum(c(0, y1))
  cumy2 <- cumsum(c(0, y2))
  cumyS <- cumsum(c(0, y1^2 + y2^2))

  cp <- rep(0, n + 1) #cp vector cp[i] = index of the last change-point for data y(1) to y(i-1)
  costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta
  index <- 0 #best index (to be found) for last change-point

  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- 1
  ###
  ### update rule Dynamic programming
  ###
  for(t in 1:n)
  {
    min_temp <- Inf
    for(i in indexSet)
    {
      val <- cumyS[t+1]-cumyS[i] - ((cumy1[t+1]-cumy1[i])^2+(cumy2[t+1]-cumy2[i])^2)/(t-i+1) + costQ[i] + beta
      if(val < min_temp)
      {
        min_temp <- val
        index <- i
      }
      costQ[t+1] <- min_temp
      cp[t+1] <- index - 1
    }

    ### PRUNING STEP
    nonpruned <- NULL
    for(i in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      if(costQ[t+1] > costQ[i] + cumyS[t+1]-cumyS[i] - ((cumy1[t+1]-cumy1[i])^2+(cumy2[t+1]-cumy2[i])^2)/(t-i+1))
      {
        nonpruned <- c(nonpruned, i)
      }
    }
    indexSet <- c(nonpruned, t+1) #add new test point
    nb[t] <- length(indexSet) ### count number of elements
  }

  cp <- cp[-1] # remove first value

  ###
  ### backtracking step
  ###
  changepoints <- n # vector of change-point to build
  current <- n

  while(current > 1)
  {
    pointval <- cp[current] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  #return(list(cp = cp, costQ = costQ, changepoints = changepoints[-1]))
  return(list(changepoints = changepoints[-1], nb = nb))
}


