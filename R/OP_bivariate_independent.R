


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
OP_2D <- function(data, beta = 4 * log(nrow(data)))
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
    for(k in 1:t)
    {
      eval <- cumyS[t+1] - cumyS[k] - ((cumy1[t+1] - cumy1[k])^2 + (cumy2[t+1] - cumy2[k])^2)/(t-k+1) + costQ[k] + beta
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
OP_2D_PELT <- function(data, beta = 4 * log(nrow(data)))
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

  cp <- rep(0, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k-1)
  costQ <- rep(0, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta

  nb <- rep(0, n) #nb element for minimization in DP
  indexSet <- 1
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:n)
  {
    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- cumyS[t+1] - cumyS[k] - ((cumy1[t+1] - cumy1[k])^2 + (cumy2[t+1] - cumy2[k])^2)/(t-k+1) + costQ[k] + beta
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[t+1] <- min_temp
    cp[t+1] <- index - 1

    ### PRUNING STEP
    nonpruned <- NULL
    for(k in indexSet)
    {
      ### PRUNING STEP: reduce indexSet using the pruning test (inquality-based)
      if(costQ[t+1] > costQ[k] + cumyS[t+1] - cumyS[k] - ((cumy1[t+1] - cumy1[k])^2 + (cumy2[t+1] - cumy2[k])^2)/(t-k+1))
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
  #return(list(cp = cp, costQ = costQ, changepoints = changepoints[-1]))
  #return(list(changepoints = changepoints[-1], nb = nb))
  return(list(changepoints = changepoints[-1], nb = nb, costQ = costQ))
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
OP_2D_1C <- function(data, beta = 4 * log(nrow(data)))
{
  ###
  ### INNER FUNCTIONS (EVAL)
  ###
  #########
  eval_var <- function(j, k)
  {
    if(j + 1 == k){return(0)}
    len <- k - j
    eval <- (cumyS[k] - cumyS[j])/len - ((cumy1[k] - cumy1[j])^2 + (cumy2[k] - cumy2[j])^2)/(len^2)
    return(eval)
  }
  #########
  eval_q_min <- function(k, t)
  {
    if(k + 1 == t){return(costQ[k] + beta)}
    eval <- cumyS[t] - cumyS[k] - ((cumy1[t] - cumy1[k])^2 + (cumy2[t] - cumy2[k])^2)/(t-k) + costQ[k] + beta
    return(eval)
  }
  #########
  eval_q <- function(k, t, t1, t2)
  {
    y1 <- cumy1[t]-cumy1[k]
    y2 <- cumy2[t]-cumy2[k]
    len <- t - k
    eval <- len*((t1 - (y1/len))^2 +  (t2 - (y2/len))^2) + cumyS[t] - cumyS[k] - ((y1)^2 + (y2)^2)/len + costQ[k] + beta
    return(eval)
  }
  #########
  eval_q_intersection <- function(j, k, t) #j<k<t
  {
    R <- (costQ[k]-costQ[j])/(k-j) - eval_var(j,k)
    t1jt <- (cumy1[t]-cumy1[k])/(t-k)
    t2jt <- (cumy2[t]-cumy2[k])/(t-k)
    t1jk <- (cumy1[k]-cumy1[j])/(k-j)
    t2jk <- (cumy2[k]-cumy2[j])/(k-j)
    D <- (t1jt - t1jk)^2 + (t2jt - t2jk)^2
    eval <- eval_q_min(k, t) + (t-k)*(sqrt(D)-sqrt(R))^2
    return(eval)
  }
  #########
  ###
  ### INITIALIZATION
  ###
  y1 <- data$y1
  y2 <- data$y2
  n <- length(y1)
  cumy1 <- cumsum(c(0, y1))
  cumy2 <- cumsum(c(0, y2))
  cumyS <- cumsum(c(0, y1^2 + y2^2))

  cp <- rep(0, n) #cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta #costQ[2] = 0

  nb <- rep(0, n) #nb element for minimization in DP
  nrows <- rep(0, n) #nb rows in info dataframe

  info <- data.frame(matrix(ncol = 3, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "m")

  ###result for first datapoint
  indexSet <- 1 #start with q1
  info[1,] <-  c(1, 1, 0)

  ###
  ### update rule Dynamic Programming
  ###
  for(t in 2:n)
  {
    ######
    ###### STEP pruning info by costQ[t] + beta = m_{t-1} + beta (costQ[1] = m_0 = -beta)
    ######
    ######### 1 ######### find omega_t^k
    omega_t_k <- rep(-Inf, length(indexSet))
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
        t1k <- (cumy1[t]-cumy1[k])/(t-k)
        t2k <- (cumy2[t]-cumy2[k])/(t-k)
        ind <- which(indexSet == k)
        if(eval_q(k, t, t1k, t2k) > eval_q(j, t, t1k, t2k)){omega_t_k[ind] <- max(omega_t_k[ind], info$m[i])}
        t1j <- (cumy1[t]-cumy1[j])/(t-j)
        t2j <- (cumy2[t]-cumy2[j])/(t-j)
        ind <- which(indexSet == j)
        if(eval_q(j, t, t1j, t2j) > eval_q(k, t, t1j, t2j)){omega_t_k[ind] <- max(omega_t_k[ind], info$m[i])}
        #remove or not?
      }
    }
    ######### 2 ######### update indexSet (removing pruned indices)
    nonpruned <- NULL
    for(k in 1:length(indexSet))
    {
      if(omega_t_k[k] < costQ[t] + beta){nonpruned <- c(nonpruned, indexSet[k])}
    }
    indexSet <- nonpruned

    ######### 3 ######### update info (removing rows with pruned indices)
    info <- info[(info$k %in% indexSet) & (info$j %in% indexSet), ]

    nb[t] <- length(indexSet) ### count number of elements for DP
    nrows[t] <- nrow(info) ### count number of rowd in info dataframe

    ###
    ### STEP updating info with new data point y(t)
    ###
    ######### 1 ######### updating 3-point already in info
    for(i in 1:nrow(info)) # update m
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
    ######### 2 ######### adding new 3-points with t

    info <- rbind(info, c(t, t, eval_q_min(t, t+1)))
    for(k in indexSet)
    {
      info <- rbind(info, c(t, k, eval_q_intersection(k, t, t+1)))
    }
    indexSet <- c(indexSet, t)


    ######
    ###### STEP find the argmin, the last change-point
    ######
    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- eval_q_min(k, t+1)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[t+1] <- min_temp
    cp[t] <- index - 1
  }

  ######
  ###### backtracking step
  ######
  changepoints <- n # vector of change-point to build
  current <- n

  while(current > 1)
  {
    pointval <- cp[current] #new last change
    changepoints <- c(pointval, changepoints) # update vector
    current <- pointval
  }
  return(list(changepoints = changepoints[-1], nb = nb, cp = cp, costQ = costQ, nrows = nrows))
}



