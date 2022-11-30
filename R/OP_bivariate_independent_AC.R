
#############################################################################################################################################
#############################################################################################################################################

###############################################
#############     OP_2D_AC    #################
###############################################

#' OP_2D_AC
#' @description Optimal Partitioning algorithm for bivariate independent time series (with PAC algorithm)
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
#' @examples
#' OP_2D_AC(dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4)))
OP_2D_AC <- function(data, beta = 4 * log(nrow(data)))
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
  eval_q_min <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(k == t){return(costQ[shift(k-1)] + beta)} ###costQ[shift(k)-1] = m_{k-1}
    return((t-k+1)*eval_var(k,t) + costQ[shift(k-1)] + beta)
  }

  #########
  eval_var <- function(j, k) ###Variance for datapoint y_j to y_{k}
  {
    if(j == k){return(0)}
    return(eval_meanS(j,k) - (eval_meany1(j,k)^2 + eval_meany2(j,k)^2))
  }
  #########
  eval_q_point <- function(k, t, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
  {
    return((t - k + 1)*((t1 - eval_meany1(k,t))^2 +  (t2 - eval_meany2(k,t))^2) + eval_q_min(k,t))
  }
  #########
  eval_q_0 <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(k == t){return(list(p = c(eval_meany1(k,t),eval_meany2(k,t)),
                           m = costQ[shift(k-1)] + beta))} ###costQ[shift(k)-1] = m_{k-1}
    return(list(p = c(eval_meany1(k,t),eval_meany2(k,t)),
                m = (t - k + 1)*eval_var(k,t) + costQ[shift(k-1)] + beta))
  }
  #########
  eval_q_1 <- function(j, k, t) ###value of m_{t}^{jk}
  {
    R <- sqrt((costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval_var(j,k-1))
    t1kt <- eval_meany1(k,t)
    t2kt <- eval_meany2(k,t)
    t1jk <- eval_meany1(j,k-1)
    t2jk <- eval_meany2(j,k-1)
    D <- sqrt((t1kt - t1jk)^2 + (t2kt - t2jk)^2)
    if(D == 0){return(list(p = c(t1kt + R, t2kt),
                           m = eval_q_min(k, t) + (t - k + 1)*R^2))}
    s <- R/D
    return(list(p = c(s*t1kt + (1-s)*t1jk, s*t2kt + (1-s)*t2jk),
                m =  eval_q_min(k, t) + (t - k + 1)*(D - R)^2))
  }
  #########
  eval_q_21 <- function(i, j, k, t) ###value of m_{t}^{ijk}
  {
    Ri2 <- (costQ[shift(k-1)] - costQ[shift(i-1)])/(k-i) - eval_var(i,k-1)
    Rj2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval_var(j,k-1)
    if((Ri2 < 0) || (Rj2 < 0)){return(list(m = Inf, p = c(0,0)))}
    Ri <- sqrt(Ri2)
    Rj <- sqrt(Rj2)
    t1i <- eval_meany1(i,k-1)
    t2i <- eval_meany2(i,k-1)
    t1j <- eval_meany1(j,k-1)
    t2j <- eval_meany2(j,k-1)
    D <- sqrt((t1i - t1j)^2 + (t2i - t2j)^2)
    if((D > Ri + Rj) | (D < Ri - Rj) | (D < Rj - Ri)){return(list(m = Inf, p = c(0,0)))}
    rho <- (1/2)*(Ri^2 - Rj^2)/D^2
    mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))
    t1A <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) + mu*(t2i - t2j)
    t2A <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) - mu*(t1i - t1j)
    mA <- eval_q_point(k, t, t1A, t2A)
    return(list(p = c(t1A, t2A), m = mA))
  }

  #########
  eval_q_22 <- function(i, j, k, t) ###value of m_{t}^{ijk}
  {
    Ri2 <- (costQ[shift(k-1)] - costQ[shift(i-1)])/(k-i) - eval_var(i,k-1)
    Rj2 <- (costQ[shift(k-1)] - costQ[shift(j-1)])/(k-j) - eval_var(j,k-1)
    if((Ri2 < 0) || (Rj2 < 0)){return(list(m = Inf, p = c(0,0)))}
    Ri <- sqrt(Ri2)
    Rj <- sqrt(Rj2)
    t1i <- eval_meany1(i,k-1)
    t2i <- eval_meany2(i,k-1)
    t1j <- eval_meany1(j,k-1)
    t2j <- eval_meany2(j,k-1)
    D <- sqrt((t1i - t1j)^2 + (t2i - t2j)^2)
    if((D > Ri + Rj) | (D < Ri - Rj) | (D < Rj - Ri)){return(list(m = Inf, p = c(0,0)))}
    rho <- (1/2)*(Ri^2 - Rj^2)/D^2
    mu <- (1/(2*D^2))*sqrt(((Ri + Rj)^2 - D^2)*(D^2 - (Ri - Rj)^2))
    t1B <- (1/2)*(t1i + t1j) - rho*(t1i - t1j) - mu*(t2i - t2j)
    t2B <- (1/2)*(t2i + t2j) - rho*(t2i - t2j) + mu*(t1i - t1j)
    mB <- eval_q_point(k, t, t1B, t2B)
    return(list(p = c(t1B, t2B), m = mB))
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

  info <- data.frame(matrix(ncol = 7, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "i", "t1", "t2", "m", "seen")
  info[1,] <-  c(1,1,1,eval_meany1(1,1),eval_meany2(1,1),0,1)

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

    # use dplyr !!!
    omega_t_k <- rep(Inf, length(indexSet)) #indexSet[u] = l, omega_t_k[u] for k = l

    for(l in 1:nrow(info))
    {
      k <- info$k[l]
      j <- info$j[l]
      i <- info$i[l]
      if(info$seen[l] > 0)
      {
        ind <- which(indexSet == k)
        omega_t_k[ind] <- min(omega_t_k[ind], info$m[l])
        ind <- which(indexSet == j)
        omega_t_k[ind] <- min(omega_t_k[ind], info$m[l])
        ind <- which(indexSet == i)
        omega_t_k[ind] <- min(omega_t_k[ind], info$m[l])
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

    test <- (info$k %in% indexSet) & (info$j %in% indexSet) & (info$i %in% indexSet)
    info <- info[test, ]

    nb[t] <- length(indexSet) ### count number of elements for DP
    nrows[t] <- nrow(info) ### count number of rows in info dataframe

    #########
    ######### STEP updating info with new data point y_{t+1}
    #########

    ######### 1 #########
    ######### 1 ######### updating 3-point already in info
    ######### 1 #########

    for(l in 1:nrow(info)) # update m with new index (time step) t+1
    {
      k <- info$k[l]
      j <- info$j[l]
      i <- info$i[l]
      if((i == j) & (j == k))
      {
        update <- eval_q_0(k, t+1)
        info$t1[l] <- update$p[1]
        info$t2[l] <- update$p[2]
        info$m[l] <- update$m
      }
      if((i == j) & (j != k))
      {
        update <- eval_q_1(j, k, t+1)
        info$t1[l] <- update$p[1]
        info$t2[l] <- update$p[2]
        info$m[l] <- update$m
      }
      if((i != j) & (j != k))
      {
        if(info$seen[l] == 1)
        {
          update <- eval_q_21(i, j, k, t+1)
          info$t1[l] <- update$p[1]
          info$t2[l] <- update$p[2]
          info$m[l] <- update$m
        }
        if(info$seen[l] == 2)
        {
          update <- eval_q_22(i, j, k, t+1)
          info$t1[l] <- update$p[1]
          info$t2[l] <- update$p[2]
          info$m[l] <- update$m
        }
      }
    }



    ######### 2 #########
    ######### 2 ######### adding new 3-points with t
    ######### 2 #########

    ##### ##### ##### ORDER 0 ##### ##### #####
    ##### ##### ##### ORDER 0 ##### ##### #####
    update <- eval_q_0(t+1, t+1)
    info <- rbind(info, c(t+1, t+1, t+1,
                          update$p[1], update$p[2],
                          update$m, 1)) #min of q_{t+1}^{t+1}

    ##### ##### ##### ORDER 1 ##### ##### #####
    ##### ##### ##### ORDER 1 ##### ##### #####

    for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint
    {
      update <- eval_q_1(j, t+1, t+1)
      info <- rbind(info, c(t+1, j, j,
                            update$p[1], update$p[2],
                            update$m, 1))
    }

    ##### ##### ##### ORDER 2 ##### ##### #####
    ##### ##### ##### ORDER 2 ##### ##### #####


    for(i in indexSet) #m_{t+1}^(ij(t+1)) optimization under two constraints
    {
      for(j in indexSet)
      {
        if(i < j)
        {
          update <- eval_q_21(i, j, t+1, t+1)
          if(update$m != Inf)
          {
            info <- rbind(info, c(t+1, j, i,
                                  update$p[1], update$p[2],
                                  update$m, 1))
          }
          if(update$m != Inf)
          {
            update <- eval_q_22(i, j, t+1, t+1)
            info <- rbind(info, c(t+1, j, i,
                                  update$p[1], update$p[2],
                                  update$m, 2))
          }
        }
      }
    }

    ###### TO DO PRUNING seen

    info$seen[info$seen <= 1] <- 1
    for(l in 1:nrow(info)) # update m with new index (time step) t+1
    {
      for(f in indexSet)
      {
        if((f != info$k[l]) & (f != info$j[l]) & (f != info$i[l]))
        {
          if(eval_q_point(f,t+1,info$t1[l],info$t2[l]) < info$m[l]){info$seen[l] <- 0}
        }
      }
    }

    test <- (info$k == info$j) | (info$j == info$i) | (info$seen > 0) #2 intersection
    info <- info[test,]
    #print("infoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfoinfo")
    #print(info)
    #print("pruning_double_point")
    #print(nrow(info))
    #pruning one constraint points
    triplePoint <- (info$k != info$j) & (info$j != info$i)
    doublePointUnseen <- (info$k != info$j) & (info$j == info$i) & (info$seen == 0)

    Index_triplePoint <- info[triplePoint,1:3]
    Index_doublePointUnSeen <- info[doublePointUnseen,1:2]

    if(nrow(Index_doublePointUnSeen) > 0){
      res <- rep(0, nrow(Index_doublePointUnSeen))

      for(f in 1:nrow(Index_doublePointUnSeen))
      {
        test12 <- any((Index_triplePoint[,1] == Index_doublePointUnSeen[f,1]) & (Index_triplePoint[,2] == Index_doublePointUnSeen[f,2]))
        test23 <- any((Index_triplePoint[,2] == Index_doublePointUnSeen[f,1]) & (Index_triplePoint[,3] == Index_doublePointUnSeen[f,2]))
        test13 <- any((Index_triplePoint[,1] == Index_doublePointUnSeen[f,1]) & (Index_triplePoint[,3] == Index_doublePointUnSeen[f,2]))
        if(test12 | test23 | test13){res[f] <- 1}
      }
      toRemove <- Index_doublePointUnSeen[as.logical(1-res),]
      #print("toRemove")
      #print(toRemove)
      if(nrow(toRemove) > 0){
        toDelete <- rep(FALSE, nrow(info))
        for(f in 1:nrow(toRemove))
        {
          toDelete <- toDelete | ((info$k == toRemove[f,1]) & (info$j == toRemove[f,2]) & (info$seen == 0))
        }
        info <- info[!toDelete,]
      }
    }
    #print(nrow(info))
    #print(info)
    #print("infoEXITinfoEXITinfoEXITinfoEXITinfoEXITinfoEXITinfoEXITinfoEXIT")
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


