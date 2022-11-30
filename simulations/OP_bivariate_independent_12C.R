
#############################################################################################################################################
#############################################################################################################################################


###############################################
#############     OP_2D_1C    #################
###############################################

#' OP_2D_1C
#' @description Optimal Partitioning algorithm for bivariate independent time series (with OP1C algorithm)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @param testMode inner pruning test modes (0 or 1)
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_1C(dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_1C <- function(data, beta = 4 * log(nrow(data)), testMode = 0)
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
    if(k == t){return(costQ[shift(k-1)] + beta)} ###costQ[shift(k)-1] = m_{k-1}
    return((t - k + 1)*eval_var(k,t) + costQ[shift(k-1)] + beta)
  }
  #########
  eval_q <- function(k, t, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
  {
    return((t - k + 1)*((t1 - eval_meany1(k,t))^2 +  (t2 - eval_meany2(k,t))^2) + eval_q_min2(k,t))
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
    return(eval_q_min2(k, t) + (t - k + 1)*(D - R)^2)
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
        if(eval_q_min2(k, t) > eval_q(j, t, t1k, t2k)){omega_t_k[ind_k] <- max(omega_t_k[ind_k], info$m[i])}
        t1j <- eval_meany1(j,t)
        t2j <- eval_meany1(j,t)
        ind_j <- which(indexSet == j)
        if(eval_q_min2(j, t) > eval_q(k, t, t1j, t2j)){omega_t_k[ind_j] <- max(omega_t_k[ind_j], info$m[i])}
      }
    }

    ######### 2 #########
    ######### 2 ######### pruning indexSet (removing pruned indices)
    ######### 2 #########

    nonpruned <- NULL
    for(k in 1:length(indexSet))
    {
      if(omega_t_k[k] <= costQ[shift(t)] + beta){nonpruned <- c(nonpruned, indexSet[k])}
    }
    indexSet <- nonpruned

    ######### 3 #########
    ######### 3 ######### pruning info (removing rows with pruned indices)
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
        info$m[i] <- eval_q_min2(k, t+1)
      }
      else
      {
        info$m[i] <- eval_q_intersection(j, k, t+1)
      }
    }

    ######### 2 #########
    ######### 2 ######### adding new 3-points with t
    ######### 2 #########

    info <- rbind(info, c(t+1, t+1, eval_q_min2(t+1, t+1))) #min of q_{t+1}^{t+1}

    ######### ######### ######### TESTMODE
    ######### ######### ######### TESTMODE
    if(testMode == 1) ### DANGER: double loop on indexSet
    {
      for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint solution
      {
        inclus <- FALSE
        for(k in indexSet)
        {
          if(j != k){if(test_inclusion(j, k, t+1)){inclus <- TRUE}}
        }
        if(inclus == FALSE){info <- rbind(info, c(t+1, j, eval_q_intersection(j, t+1, t+1)))}
      }
    }
    if(testMode == 0)
    {
      for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint solution
      {
        info <- rbind(info, c(t+1, j, eval_q_intersection(j, t+1, t+1)))
      }
    }
    ######### ######### ######### TESTMODE
    ######### ######### ######### TESTMODE

    ###add new index t
    indexSet <- c(indexSet, t+1)

    #########
    ######### STEP find the argmin, the last change-point
    #########

    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- eval_q_min2(k, t+1) ### find min of q_{t+1}^k
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


#############################################################################################################################################
#############################################################################################################################################


###############################################
#############     OP_2D_2C    #################
###############################################

#' OP_2D_2C
#' @description Optimal Partitioning algorithm for bivariate independent time series (with OP2C algorithm)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @param testMode inner pruning test modes (0, 1 or 2)
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_2C(dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_2C <- function(data, beta = 4 * log(nrow(data)), testMode = 2)
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
    if(k == t){return(costQ[shift(k-1)] + beta)} ###costQ[shift(k)-1] = m_{k-1}
    return((t - k + 1)*eval_var(k,t) + costQ[shift(k-1)] + beta)
  }
  #########
  eval_q <- function(k, t, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
  {
    return((t - k + 1)*((t1 - eval_meany1(k,t))^2 +  (t2 - eval_meany2(k,t))^2) + eval_q_min2(k,t))
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
    return(eval_q_min2(k, t) + (t - k + 1)*(D - R)^2)
  }
  #########
  position <- function(j, k, t)
  {
    R <- sqrt((costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval_var(j,k-1))
    t1kt <- eval_meany1(k,t)
    t2kt <- eval_meany2(k,t)
    t1jk <- eval_meany1(j,k-1)
    t2jk <- eval_meany2(j,k-1)
    D <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)
    if(D == 0){return(c(t1kt + R, t2kt))}
    s <- R/D
    return(c(s*t1kt + (1-s)*t1jk, s*t2kt + (1-s)*t2jk))
  }
  #########
  eval_q_2intersection <- function(i, j, k, t) ###value of m_{t}^{ijk}
  {
    Ri <- sqrt((costQ[shift(k-1)] - costQ[shift(i-1)])/(k-i) - eval_var(i,k-1))
    Rj <- sqrt((costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval_var(j,k-1))
    t1i <- eval_meany1(i,k-1)
    t2i <- eval_meany2(i,k-1)
    t1j <- eval_meany1(j,k-1)
    t2j <- eval_meany2(j,k-1)
    D <- sqrt((t1i - t1j)^2 + (t2i - t2j)^2)
    if((D > Ri + Rj) | (D < Ri - Rj) | (D < Rj - Ri)){return(list(m = -Inf, p = c(0,0)))}
    rho <- (1/2)*(Ri^2 - Rj^2)/D^2
    mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))
    t1A <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) + mu*(t2i - t2j)
    t2A <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) - mu*(t1i - t1j)
    t1B <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) - mu*(t2i - t2j)
    t2B <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) + mu*(t1i - t1j)
    mA <- eval_q(k, t, t1A, t2A)
    mB <- eval_q(k, t, t1B, t2B)
    res <- which.min(c(mA, mB))
    if(res == 1){return(list(m = mA, p = c(t1A, t2A)))}
    if(res == 2){return(list(m = mB, p = c(t1B, t2B)))}
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

  #########
  ###
  ### INITIALIZATION
  ###
  costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta #costQ[2] = 0
  cp <- rep(0, n) #cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  nb <- rep(0, n) #nb element for minimization in DP
  nrows <- rep(0, n) #nb rows in info dataframe
  nrows3 <- rep(0, n) #nb rows in info dataframe

  info <- data.frame(matrix(ncol = 3, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "m")
  info[1,] <-  c(1, 1, 0)

  info3 <- data.frame(matrix(ncol = 10, nrow = 0)) ### info for pruning, triple points (m_t^{ijk})
  colnames(info3) <- c("k", "j", "i", "ti1", "ti2", "tj1", "tj2", "ijk1", "ijk2", "m")
  info3[1,] <-  rep(0,10)

  ###result for first datapoint
  indexSet <- 1 #start with q1

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
        if(eval_q_min2(k, t) > eval_q(j, t, t1k, t2k)){omega_t_k[ind_k] <- max(omega_t_k[ind_k], info$m[i])}
        t1j <- eval_meany1(j,t)
        t2j <- eval_meany1(j,t)
        ind_j <- which(indexSet == j)
        if(eval_q_min2(j, t) > eval_q(k, t, t1j, t2j)){omega_t_k[ind_j] <- max(omega_t_k[ind_j], info$m[i])}
      }
    }

    ################# TRIPLE INTERSECTION
    if(nrow(info3) > 1)
    {
      for(l in 2:nrow(info3)) #remove row 1 = 0...0
      {
        k <- info3$k[l]
        j <- info3$j[l]
        i <- info3$i[l]
        ti1 <- info3$ti1[l]
        ti2 <- info3$ti2[l]
        tj1 <- info3$tj1[l]
        tj2 <- info3$tj2[l]
        ind_k <- which(indexSet == k)
        test <- FALSE
        if((eval_q(k, t, ti1, ti2) > eval_q(j, t, ti1, ti2)) && (eval_q(k, t, tj1, tj2) > eval_q(i, t, tj1, tj2))){test <- TRUE}
        if(test == TRUE){omega_t_k[ind_k] <- max(omega_t_k[ind_k], info3$m[l])}
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
    ######### 3 ######### pruning info (removing rows with pruned indices)
    ######### 3 #########

    info <- info[(info$k %in% indexSet) & (info$j %in% indexSet), ]
    test <- (info3$k %in% indexSet) & (info3$j %in% indexSet) & (info3$i %in% indexSet)
    test[1] <- TRUE #save the first row
    info3 <- info3[test, ]

    ######### ######### ######### TESTMODE
    ######### ######### ######### TESTMODE
    if(testMode == 2)
    {
      if(nrow(info3) > 1)
      {
        stay <- rep(TRUE, nrow(info3))
        for(l in 2:nrow(info3))
        {
          for(m in indexSet)
          {
            if((m != info3$k[l]) && (m != info3$j[l]) && (m != info3$i[l]))
            {if(eval_q(m, t, info3$ijk1[l], info3$ijk2[l]) < info3$m[l]){stay[l] <- FALSE}}
          }
        }
        info3 <- info3[stay, ]
      }
    }
    ######### ######### ######### TESTMODE
    ######### ######### ######### TESTMODE

    ### + IF index not in info3, remove (with pb for 2 quadratics) from info

    nb[t] <- length(indexSet) ### count number of elements for DP
    nrows[t] <- nrow(info) ### count number of rows in info dataframe
    nrows3[t] <- nrow(info3) ### count number of rows in info3 dataframe

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
        info$m[i] <- eval_q_min2(k, t+1)
      }
      else
      {
        info$m[i] <- eval_q_intersection(j, k, t+1)
      }
    }


    ################# TRIPLE INTERSECTION OK
    ################# TRIPLE INTERSECTION OK
    if(nrow(info3) > 1)
    {
      for(l in 2:nrow(info3))
      {
        k <- info3$k[l]
        j <- info3$j[l]
        i <- info3$i[l]
        position_i <- position(i,k,t+1)
        position_j <- position(j,k,t+1)
        info3$ti1[l] <- position_i[1]
        info3$ti2[l] <- position_i[2]
        info3$tj1[l] <- position_j[1]
        info3$tj2[l] <- position_j[2]
        res <- eval_q_2intersection(i,j,k,t+1)
        info3$ijk1[l] <- res$p[1]
        info3$ijk2[l] <- res$p[2]
        info3$m[l] <- res$m
      }
    }


    ######### 2 #########
    ######### 2 ######### adding new 3-points with t
    ######### 2 #########

    info <- rbind(info, c(t+1, t+1, eval_q_min2(t+1, t+1))) #min of q_{t+1}^{t+1}

    ######### ######### ######### TESTMODE
    ######### ######### ######### TESTMODE
    if((testMode == 1) || (testMode == 2)) ### DANGER: double loop on indexSet
    {
      for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint solution
      {
        inclus <- FALSE
        for(k in indexSet)
        {
          if(j != k){if(test_inclusion(j, k, t+1)){inclus <- TRUE}}
        }
        if(inclus == FALSE){info <- rbind(info, c(t+1, j, eval_q_intersection(j, t+1, t+1)))}
      }
    }
    if(testMode == 0)
    {
      for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint solution
      {
        info <- rbind(info, c(t+1, j, eval_q_intersection(j, t+1, t+1)))
      }
    }
    ######### ######### ######### TESTMODE
    ######### ######### ######### TESTMODE


    ################# TRIPLE INTERSECTION OK
    ################# TRIPLE INTERSECTION OK
    for(i in indexSet)
    {
      for(j in indexSet)
      {
        if(i< j)
        {
          res <- eval_q_2intersection(i,j,t+1,t+1)
          info3 <- rbind(info3,
                         c(t+1, j, i,
                           position(i,t+1,t+1),
                           position(j,t+1,t+1),
                           res$p, res$m))
        }
      }
    }

    ###add new index t
    indexSet <- c(indexSet, t+1)

    #########
    ######### STEP find the argmin, the last change-point
    #########

    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- eval_q_min2(k, t+1) ### find min of q_{t+1}^k
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
  return(list(changepoints = changepoints[-1],  nb = nb[1:(n-1)], nrows = nrows[1:(n-1)], nrows3 = nrows3[1:(n-1)]))
}



