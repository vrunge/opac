
############################################
#############     OP_Reg    ################
############################################

#' OP_Reg
#' @description Optimal Partitioning algorithm for changes in simple regression (no pruning)
#' @param data a dataframe with two components: x and y, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment)
#' @examples
#' OP_Reg(dataGeneratorRegression(chpts = c(50,100,150), A = c(-1,1,-1), B = c(-1,1,-1)))
OP_Reg <- function(data, beta = 4 * log(nrow(data)))
{
  ###
  ### INITIALIZATION
  ###
  x <- data$x
  y <- data$y
  n <- length(y)
  cumX <- cumsum(c(0, x))
  cumY <- cumsum(c(0, y))
  cumXY <- cumsum(c(0, x*y))
  cumSX <- cumsum(c(0, x^2))
  cumSY <- cumsum(c(0, y^2))

  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  index <- 0 #best index (to be found) for last change-point

  ###
  ### UPDATE rule Dynamic Programming
  ###
  for(t in 1:n)
  {
    min_temp <- Inf
    index <- 1
    for(k in 1:t)
    {
      if(t == k){eval <- Inf}
      else
      {
        t1 <- ((t-k+1)*(cumXY[t+1]-cumXY[k]) - (cumX[t+1]-cumX[k])*(cumY[t+1]-cumY[k]))/((t-k+1)*(cumSX[t+1]-cumSX[k]) - (cumX[t+1]-cumX[k])^2)
        t2 <- ((cumSX[t+1]-cumSX[k])*(cumY[t+1]-cumY[k]) - (cumX[t+1]-cumX[k])*(cumXY[t+1]-cumXY[k]))/((t-k+1)*(cumSX[t+1]-cumSX[k]) - (cumX[t+1]-cumX[k])^2)
        eval <-  (cumSX[t+1]-cumSX[k])*t1^2 + 2*(cumX[t+1]-cumX[k])*t1*t2 + (t-k+1)*t2^2
        eval <- eval - 2*(cumY[t+1]-cumY[k])*t2 - 2*(cumXY[t+1]-cumXY[k])*t1 + cumSY[t+1] - cumSY[k] + costQ[k] + beta
      }
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[t+1] <- min_temp
    cp[t+1] <- index - 1
  }

  ###
  ### backtracking step
  ###
  cp <- cp[-1] # remove first value
  changepoints <- n # vector of change-point to build
  current <- n
  print(cp)
  while(current > 1)
  {
    pointval <- cp[current] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  #return(list(changepoints = changepoints[-1]))
  return(list(cp = cp, costQ = costQ, changepoints = changepoints[-1]))
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
#' OP_Reg_PELT(dataGeneratorRegression(chpts = c(50,100,150), A = c(-1,1,-1), B = c(-1,1,-1)))
OP_Reg_PELT <- function(data, beta = 4 * log(nrow(data)))
{
  ###
  ### INNER FUNCTIONS (EVAL)
  ###
  #########
  eval_segment_kt <- function(k, t)
  {
    eval <- 0
    if(k < t)
    {
      t1 <- ((t-k+1)*(cumXY[t+1]-cumXY[k]) - (cumX[t+1]-cumX[k])*(cumY[t+1]-cumY[k]))/((t-k+1)*(cumSX[t+1]-cumSX[k]) - (cumX[t+1]-cumX[k])^2)
      t2 <- ((cumSX[t+1]-cumSX[k])*(cumY[t+1]-cumY[k]) - (cumX[t+1]-cumX[k])*(cumXY[t+1]-cumXY[k]))/((t-k+1)*(cumSX[t+1]-cumSX[k]) - (cumX[t+1]-cumX[k])^2)
      eval <-  (cumSX[t+1]-cumSX[k])*t1^2 + 2*(cumX[t+1]-cumX[k])*t1*t2 + (t-k+1)*t2^2
      eval <- eval - 2*(cumY[t+1]-cumY[k])*t2 - 2*(cumXY[t+1]-cumXY[k])*t1 + cumSY[t+1] - cumSY[k]
    }
    return(eval)
  }

  ###
  ### INITIALIZATION
  ###
  x <- data$x
  y <- data$y
  n <- length(y)
  cumX <- cumsum(c(0, x))
  cumY <- cumsum(c(0, y))
  cumXY <- cumsum(c(0, x*y))
  cumSX <- cumsum(c(0, x^2))
  cumSY <- cumsum(c(0, y^2))

  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta

  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- 1
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
    costQ[t+1] <- min_temp
    cp[t+1] <- index - 1

    ### PRUNING STEP
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      if(costQ[t+1] > costQ[k] + eval_segment_kt(k,t))
      {
        nonpruned <- c(nonpruned, k)
      }
    }
    indexSet <- c(nonpruned, t+1) #add new test point
    nb[t] <- length(indexSet) ### count number of elements
  }

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
  #return(list(changepoints = changepoints[-1]))
  return(list(changepoints = changepoints[-1], nb = nb[1:(n-1)]))
}




