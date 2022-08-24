

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

#############################################################################################################################################
#############################################################################################################################################


###############################################
#############     OP_2D_1C    #################
###############################################

#' OP_2D_1C
#' @description Optimal Partitioning algorithm for bivariate independent time series (with OP1C algorithm)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_1C(dataGenerator2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_1C <- function(data, beta = 4 * log(nrow(data)))
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
    if(k == t){return(costQ[shift(k-1)] + beta)} ###costQ[shift(k)-1] = m_{k-1}
    return((t - k + 1)*eval_var(k,t) + costQ[shift(k-1)] + beta)
  }
  #########
  eval_q <- function(k, t, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
  {
    return((t - k + 1)*((t1 - eval_meany1(k,t))^2 +  (t2 - eval_meany2(k,t))^2) + eval_q_min(k,t))
  }
  #########
  eval_q_intersection <- function(j, k, t) ###value of m_{t}^{jk}
  {
    R <- sqrt((costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval_var(j,k-1))
    t1kt <- eval_meany1(k,t)
    t2kt <- eval_meany2(k,t)
    t1jk <- eval_meany1(j,k-1)
    t2jk <- eval_meany2(j,k-1)
    D <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)
    return(eval_q_min(k, t) + (t - k + 1)*(D - R)^2)
  }

  #########  #########  #########
  #########  #########  #########
  #########  #########  #########

  test_inclusion <- function(j, k, t) #j < k to prune j
  {
    Rsmall <- sqrt((costQ[shift(t-1)] - costQ[shift(j-1)])/(t-j) - eval_var(j,t-1))
    Rbig <- sqrt((costQ[shift(t-1)] - costQ[shift(k-1)])/(t-k) - eval_var(k,t-1))
    t1jt <- eval_meany1(j,t-1)
    t2jt <- eval_meany2(j,t-1)
    t1kt <- eval_meany1(k,t-1)
    t2kt <- eval_meany2(k,t-1)
    D <- sqrt((t1kt - t1jt)^2 + (t2kt - t2jt)^2)
    return((D + Rsmall) <= Rbig)
  }

  test_intersection <- function(j, k, t)
  {
    Rsmall <- sqrt((costQ[shift(t-1)] - costQ[shift(j-1)])/(t-j) - eval_var(j,t-1))
    Rbig <- sqrt((costQ[shift(t-1)] - costQ[shift(k-1)])/(t-k) - eval_var(k,t-1))
    t1jt <- eval_meany1(j,t-1)
    t2jt <- eval_meany2(j,t-1)
    t1kt <- eval_meany1(k,t-1)
    t2kt <- eval_meany2(k,t-1)
    D <- sqrt((t1kt - t1jt)^2 + (t2kt - t2jt)^2)
    return(D > (Rsmall + Rbig))
  }

  #########
  ###
  ### INITIALIZATION
  ###
  costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta #costQ[2] = 0
  cp <- rep(0, n) #cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  nb <- rep(0, n) #nb element for minimization in DP
  nrows <- rep(0, n) #nb rows in info dataframe

  info <- data.frame(matrix(ncol = 3, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "m")

  ###result for first datapoint
  indexSet <- 1 #start with q1
  info[1,] <-  c(1, 1, 0)

  #########
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:(n-1)) # at t, transform Q_t into Q_{t+1}
  {
    #########
    ######### STEP pruning info by costQ[shift(t)]  + beta = m_{t} + beta
    ######### = min Q_{t}(.,.) + beta
    ######### (costQ[shift(0)] = costQ[1] = m_0 = -beta)
    #########

    ######### 1 #########
    ######### 1 ######### find omega_{t}^k
    ######### 1 #########

    omega_t_k <- rep(-Inf, length(indexSet)) #indexSet[u] = l, omega_t_k[u] for k = l
    for(i in 1:nrow(info))
    {
      k <- info$k[i]
      j <- info$j[i]
      if(j == k)
      {
        ind <- which(indexSet == k)
        omega_t_k[ind] <- max(omega_t_k[ind], info$m[i])
      }
      else
      {
        t1k <- eval_meany1(k,t)
        t2k <- eval_meany2(k,t)
        ind_k <- which(indexSet == k)
        if(eval_q_min(k, t) > eval_q(j, t, t1k, t2k)){omega_t_k[ind_k] <- max(omega_t_k[ind_k], info$m[i])}
        t1j <- eval_meany1(j,t)
        t2j <- eval_meany1(j,t)
        ind_j <- which(indexSet == j)
        if(eval_q_min(j, t) > eval_q(k, t, t1j, t2j)){omega_t_k[ind_j] <- max(omega_t_k[ind_j], info$m[i])}
      }
    }

    ######### 2 #########
    ######### 2 ######### update indexSet (removing pruned indices)
    ######### 2 #########

    nonpruned <- NULL
    for(k in 1:length(indexSet))
    {
      if(omega_t_k[k] <= costQ[shift(t)] + beta){nonpruned <- c(nonpruned, indexSet[k])}
    }
    indexSet <- nonpruned

    ######### 3 #########
    ######### 3 ######### update info (removing rows with pruned indices)
    ######### 3 #########

    info <- info[(info$k %in% indexSet) & (info$j %in% indexSet), ]
    nb[t] <- length(indexSet) ### count number of elements for DP
    nrows[t] <- nrow(info) ### count number of rowd in info dataframe

    #########
    ######### STEP updating info with new data point y_{t+1}
    #########

    ######### 1 #########
    ######### 1 ######### updating 3-point already in info
    ######### 1 #########

    for(i in 1:nrow(info)) # update m with new index (time step) t+1
    {
      k <- info$k[i]
      j <- info$j[i]
      if(k == j)
      {
        info$m[i] <- eval_q_min(k, t+1)
      }
      else
      {
        info$m[i] <- eval_q_intersection(j, k, t+1)
      }
    }

    ######### 2 #########
    ######### 2 ######### adding new 3-points with t
    ######### 2 #########

    info <- rbind(info, c(t+1, t+1, eval_q_min(t+1, t+1))) #min of q_{t+1}^{t+1}

    for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint solution
    {
      inclus <- FALSE
      for(k in indexSet)
      {
        if(j != k){if(test_inclusion(j, k, t+1)){inclus <- TRUE}}
      }
      if(inclus == FALSE){info <- rbind(info, c(t+1, j, eval_q_intersection(j, t+1, t+1)))}
    }

    ###add new index t
    indexSet <- c(indexSet, t+1)

    #########
    ######### STEP find the argmin, the last change-point
    #########

    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- eval_q_min(k, t+1) ### find min of q_{t+1}^k
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t+1)] <- min_temp # this is find m_{t+1}
    cp[t+1] <- index - 1
  }

  #########
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
  return(list(changepoints = changepoints[-1],  nb = nb[1:(n-1)], nrows = nrows[1:(n-1)]))
}



###############################################
#############     OP_2D_2C    #################
###############################################

#' OP_2D_2C
#' @description Optimal Partitioning algorithm for bivariate independent time series (with OP2C algorithm)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_2C(dataGenerator2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_2C <- function(data, beta = 4 * log(nrow(data)))
{


}


