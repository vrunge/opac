
###############################################
#############     OP_2D_1C    #################
###############################################

#' OP_2D_1C
#'
#' @description Optimal Partitioning algorithm for bivariate independent time series (with OP1C algorithm)
#' @param data a data-frame with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_1C(dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
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
  ### INITIALIZATION
  ###
  costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta #costQ[2] = 0
  cp <- rep(0, n + 1) #cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  #cp[shift(1)] <- 0
  nb <- rep(0, n) #nb element for minimization in DP
  nrows <- rep(0, n) #nb rows in info data-frame

  info <- data.frame(matrix(ncol = 3, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "m")

  ###result for first data-point
  indexSet <- 1 #start with q1
  info[1,] <-  c(1, 1, 0)
  nb[1] <- 1
  nrows[1] <- 1

  #####################
  ###
  ### update rule Dynamic Programming
  ###
  for(t in 1:(n-1)) # at t, transform Q_t into Q_{t+1}
  {
    ######### AAAAAAAAAAAAAAAAAAAAAA #########
    #########
    ######### STEP pruning info by costQ[shift(t)]  + beta = m_{t} + beta
    ######### = min Q_{t}(.,.) + beta
    ######### (costQ[shift(0)] = costQ[1] = m_0 = -beta)
    #########

    ######### 1 ######### find omega_{t}^k
    omega_t_k <- omega_t_k_fct(t, info, nrows[t],
                               indexSet, nb[t],
                               costQ, cumy1, cumy2, cumyS, beta)

    ######### 2 ######### pruning indexSet (removing pruned indices)
    nonpruned <- NULL
    for(k in 1:nb[t]){if(omega_t_k[k] <= costQ[shift(t)] + beta){nonpruned <- c(nonpruned, indexSet[k])}}
    indexSet <- nonpruned

    ######### 3 ######### pruning info (removing rows with pruned indices)
    info <- info[(info$k %in% indexSet) & (info$j %in% indexSet), ]
    nrows[t+1] <- nrow(info) ### count number of rows in info data-frame (update later again)

    ######### BBBBBBBBBBBBBBBBBBBBBBBBBB #########
    #########
    ######### STEP updating info with new data point y_{t+1}
    #########

    ######### 1 ######### updating 3-point already in info

    for(i in 1:nrows[t+1]) # update m with new index (time step) t+1
    {
      k <- info$k[i]
      j <- info$j[i]
      if(j == k)
      {
        info$m[i] <- eval_q_min(costQ, cumy1, cumy2, cumyS, k, t+1, beta)
      }
      else
      {
        R2 <- sqrt((costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval_var(cumy1, cumy2, cumyS, j, k))
        info$m[i] <- eval_q_intersection(R2, costQ, cumy1, cumy2, cumyS, j, k, t+1, beta)
      }
    }

    ######### 2 ######### adding new 3-points with t

    info <- rbind(info, c(t+1, t+1, eval_q_min(costQ, cumy1, cumy2, cumyS, t+1, t+1, beta))) #min of q_{t+1}^{t+1}

    for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint solution
    {
      R2 <- (costQ[shift(t)] - costQ[shift(j-1)])/(t+1-j) - eval_var(cumy1, cumy2, cumyS, j, t+1)
      if(R2 > 0)
      {
        info <- rbind(info, c(t+1, j, eval_q_intersection(R2, costQ, cumy1, cumy2, cumyS, j, t+1, t+1, beta)))
      }
    }

    ### add new index t
    indexSet <- c(indexSet, t+1)

    nb[t+1] <- length(indexSet) ### count number of elements for DP
    nrows[t+1] <- nrow(info) ### count number of rows in info data-frame

    ######### CCCCCCCCCCCCCCCCCCCCCCCCCCCc #########
    #########
    ######### STEP find the argmin, the last change-point
    #########

    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- eval_q_min(costQ, cumy1, cumy2, cumyS, k, t+1, beta) ### find min of q_{t+1}^k
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t+1)] <- min_temp # this is m_{t+1}
    cp[shift(t+1)] <- index - 1
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

  return(list(changepoints = changepoints[-1],  nb = nb, nrows = nrows))
}
