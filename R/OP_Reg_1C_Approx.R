
###############################################
#############     OP_Reg_1C    ################
###############################################

#' OP_Reg_1C
#'
#' @description Optimal Partitioning algorithm for change in simple regression (with OP1C algorithm)
#' @param data a data-frame with two components: x and y, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
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
  ### INITIALIZATION
  ###
  cp <- rep(NA, n + 1) #cp vector cp[k] = index of the last change-point for data y(1) to y(k)
  cp[1] <- 0
  costQ <- rep(NA, n + 1) # costQ[k] optimal cost for data y(1) to y(k-1)
  costQ[1] <- -beta
  nb <- rep(NA, n) #nb element for minimization in DP
  nrows <- rep(NA, n) #nb rows in info data-frame

  info <- data.frame(matrix(ncol = 3, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "m")

  ###result for first data-point
  indexSet <- 1 #start with q1
  costQ[shift(1)] <- Inf #no change point at position 1
  cp[shift(1)] <- 0
  info[1,] <-  c(1, 1, Inf) #cost of one point = Inf for simple regression
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

    omega_t <-  omega_t_Reg_fct_1C(t, info, nrows[t],
                                   indexSet, nb[t],
                                   costQ, cumX, cumY, cumXY, cumSX, cumSY, beta)

    ######### 2 ######### pruning indexSet (removing pruned indices)
    nonpruned <- NULL
    for(k in 1:nb[t]){if(omega_t[k] <= costQ[shift(t)] + beta){nonpruned <- c(nonpruned, indexSet[k])}}
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
        info$m[i] <- evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t+1, beta)
      }
      else
      {
        info$m[i] <- evalReg_q_1_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t+1, beta)
      }
    }

    ######### 2 ######### adding new 3-points with t

    info <- rbind(info, c(t+1, t+1, -Inf)) #min of q_{t+1}^{t+1} no segment of length 1 and no pruning

    for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint solution
    {
      coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t+1, beta)
      coefft <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, t+1, t+1, beta)
      coeff1C <- coeffj - coefft

      if(isAnEllipse(coeff1C) == TRUE)
      {
        info <- rbind(info, c(t+1, j, -Inf))
      }
    }

    ### add new index t
    indexSet <- c(indexSet, t+1)

    nb[t+1] <- length(indexSet) ### count number of elements for DP
    nrows[t+1] <- nrow(info) ### count number of rows in info data-frame

    ######### CCCCCCCCCCCCCCCCCCCCCCCCCCCC #########
    #########
    ######### STEP find the argmin, the last change-point
    #########

    #print(info)

    min_temp <- Inf
    for(k in indexSet[-nb[t+1]]) # remove index t+1
    {
      eval <- evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t+1, beta)
      if(eval < min_temp){min_temp <- eval; index <- k}
    }
    costQ[shift(t+1)] <- min_temp
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


