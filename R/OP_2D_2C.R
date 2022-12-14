
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
#' OP_2D_2C(dataGenerator_2D(chpts = c(3,10,12), means1 = c(0,1,0), means2 = c(7,1,-4)))
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
  ### INITIALIZATION
  ###
  costQ <- rep(NA, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta
  cp <- rep(NA, n + 1) #cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  cp[1] <- 0
  nb <- rep(NA, n) #nb element for minimization in DP
  nrows <- rep(NA, n) #nb rows in info data-frame

  ### indices k,j,i; min values m and M, positions (p1,p2) and (q1,q2)
  info <- data.frame(matrix(ncol = 9, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "i", "m", "M","p1", "p2", "q1", "q2")


  ###result for first data-point
  indexSet <- 1 #start with q1
  costQ[shift(1)] = 0
  cp[shift(1)] <- 0
  info[1,] <-  c(1, NA, NA, 0, NA, y1[1], y2[1], Inf, Inf)
  nb[1] <- 1
  nrows[1] <- 1

  #########
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
    omega_t <- omega_t_fct_2C(t, info,
                              indexSet, nb[t],
                              costQ, cumy1, cumy2, cumyS, beta)

    ######### 2 ######### pruning indexSet (removing pruned indices)
    nonpruned <- NULL
    for(k in 1:nb[t]){if(omega_t[k] <= costQ[shift(t)] + beta){nonpruned <- c(nonpruned, indexSet[k])}}
    indexSet <- nonpruned

    ######### 3 #########
    ######### 3 ######### pruning info (removing rows with pruned indices)
    ######### 3 #########


    ######### 3 ######### pruning info (removing rows with pruned indices)

    prunedByPosition <- (info$p1 == Inf) & (info$q1 == Inf) ### unseen
    info <- info[((info$k %in% indexSet) &
                (info$j %in% c(indexSet, NA)) &
                (info$i %in% c(indexSet, NA))) &
                (!prunedByPosition)  , ]

    nrows[t+1] <- nrow(info) ### count number of rows in info data-frame (update later again)

    ######### BBBBBBBBBBBBBBBBBBBBBBBBBB #########
    #########
    ######### STEP updating info with new data point y_{t+1}
    #########

    ######### 1 ######### updating 3-point already in info

    for(l in 1:nrows[t+1]) # update m with new index (time step) t+1
    {
      k <- info$k[l]
      j <- info$j[l]
      i <- info$i[l]

      if(is.na(j))
      {
        res <- eval2D_q_0(costQ, cumy1, cumy2, cumyS, k, t+1, beta)
        info$p1[l] <- res$p[1]
        info$p2[l] <- res$p[2]
        info$m[l] <- res$m
      }
      if((!is.na(j)) & is.na(i))
      {
        res <- eval2D_q_1(costQ, cumy1, cumy2, cumyS, j, k, t+1, beta)
        info$p1[l] <- res$p[1]
        info$p2[l] <- res$p[2]
        info$m[l] <- res$m
      }
      if((!is.na(j)) & (!is.na(i)))
      {
        res1 <- eval2D_q_21(costQ, cumy1, cumy2, cumyS, i, j, k, t+1, beta)
        res2 <- eval2D_q_22(costQ, cumy1, cumy2, cumyS, i, j, k, t+1, beta)
        if(res1$p[1] != Inf) #update if point is seen
        {
          info$p1[l] <- res1$p[1]
          info$p2[l] <- res1$p[2]
          info$m[l] <- res1$m
        }
        if(res2$p[1] != Inf) #update if point is seen
        {
          info$q1[l] <- res2$p[1]
          info$q2[l] <- res2$p[2]
          info$M[l] <- res2$m
        }
      }
    }


    ######### 2 ######### adding new 3-points with t

    temp <- eval2D_q_0(costQ, cumy1, cumy2, cumyS, t+1, t+1, beta)
    info <- rbind(info, c(t+1, NA, NA,
                          temp$m, NA, temp$p[1], temp$p[2], Inf, Inf))

    for(j in indexSet) #m_{t+1}^(j(t+1)) optimization under one constraint
    {
      temp <- eval2D_q_1(costQ, cumy1, cumy2, cumyS, j, t+1, t+1, beta)
      info <- rbind(info, c(t+1, j, NA,
                            temp$m, NA, temp$p[1], temp$p[2], Inf, Inf))

      for(i in indexSet) #m_{t+1}^(ij(t+1)) optimization under two constraint
      {
        if((i< j) & eval2D_circleIntersection_ijk(costQ, cumy1, cumy2, cumyS, i, j, t+1))
        {
          temp1 <- eval2D_q_21(costQ, cumy1, cumy2, cumyS, i, j, t+1, t+1, beta)
          temp2 <- eval2D_q_22(costQ, cumy1, cumy2, cumyS, i, j, t+1, t+1, beta)

          if((temp1$p[1] != Inf) | (temp2$p[1] != Inf))
          {
            info <- rbind(info, c(t+1, j, i,
                                temp1$m, temp2$m,
                                temp1$p[1], temp1$p[2], temp2$p[1], temp2$p[2]))
          }
        }
      }
    }

    ###
    ###
    ### additional pruning of 3points with "point covering" using last quadratics t+1
    ###
    ###
    info <- points_toInf(info, costQ, cumy1, cumy2, cumyS, t+1, beta)

    ### add new index t
    indexSet <- c(indexSet, t+1)

    nb[t+1] <- length(indexSet) ### count number of elements for DP
    nrows[t+1] <- nrow(info) ### count number of rows in info data-frame

    ######### CCCCCCCCCCCCCCCCCCCCCCCCCCCC #########
    #########
    ######### STEP find the argmin, the last change-point
    #########

    min_temp <- Inf
    for(k in indexSet)
    {
      eval <- eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t+1, beta) ### find min of q_{t+1}^k
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



