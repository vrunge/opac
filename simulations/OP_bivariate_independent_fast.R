
#################################################
#############     OP_2D_fast    #################
#################################################

#' OP_2D_fast
#' @description Optimal Partitioning algorithm for bivariate independent time series (no pruning)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment)
#' @examples
#' OP_2D_fast(dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,5,0), means2 = c(7,1,-4)))
OP_2D_fast <- function(data, beta = 4 * log(nrow(data)))
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
  eval_q_min2 <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
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
  eval_min <- data.table(k = 1:n, min_k_t = rep(NA, n))
  index <- 0 #best index (to be found) for last change-point

  #########
  ###
  ### UPDATE rule Dynamic Programming
  ###
  for(t in 1:n) # at t, transform Q_{t-1} into Q_{t}
  {
    ########### TO DO ###########
    #rows 1 to t. update
    eval_min <- eval_min %>%
      mutate(min_k_t = -1)
    #rows 1 to t. find min and argmin (min_temp and index)
    min_temp <- 0
    ###########

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

