
###############################################
#############     OP_2D_AC    #################
###############################################

#' OP_2D_AC
#' @description Optimal Partitioning algorithm for bivariate independent time series (with OPAC algorithm)
#' @param data a data-frame with two components: y1 and y2, time series of same length
#' @param beta penalty value
#' @return a list with the change-point elements (each last index of each segment) and a vector nb counting the number of non-pruned elements at each iteration
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
  ### INITIALIZATION
  ###
  costQ <- rep(NA, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
  costQ[1] <- -beta
  cp <- rep(NA, n + 1) #cp vector cp[i] = index of the last change-point for data y(1) to y(i)
  cp[1] <- 0
  nb <- rep(NA, n) #nb element for minimization in DP
  nrows <- rep(NA, n) #nb rows in info data-frame

  ### indices k,j,i; min values m and M, positions (p1,p2) and (q1,q2)
  info <- data.frame(matrix(ncol = 7, nrow = 0)) ### info for pruning
  colnames(info) <- c("k", "j", "i", "m", "p1", "p2", "surface")


  ###result for first data-point
  indexSet <- 1 #start with q1
  costQ[shift(1)] = 0
  cp[shift(1)] <- 0
  info[1,] <-  c(1, NA, NA, 0, y1[1], y2[1], 1)
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

    omega_t <- rep(Inf, nb[t])
    infoSurface <- info[info$surface >= 1,]
    for(i in 1:nb[t])
    {
      l <- indexSet[i]
      ind_k <- which(infoSurface$k == l)
      ind_j <- which(infoSurface$j == l)
      ind_i <- which(infoSurface$i == l)
      omega_t[i] <- min(omega_t[i], infoSurface$m[c(ind_k, ind_j, ind_i)])
    }

    ######### 2 ######### pruning indexSet (removing pruned indices)
    nonpruned <- NULL
    for(k in 1:nb[t]){if(omega_t[k] <= costQ[shift(t)] + beta){nonpruned <- c(nonpruned, indexSet[k])}}
    indexSet <- nonpruned

    ######### 3 #########
    ######### 3 ######### pruning info (removing rows with pruned indices)
    ######### 3 #########

    ######### 3 ######### pruning info (removing rows with pruned indices)

    info <- info[(info$k %in% indexSet) &
                    (info$j %in% c(indexSet, NA)) &
                    (info$i %in% c(indexSet, NA)), ]

    nrows[t+1] <- nrow(info) ### count number of rows in info data-frame (update later again)

    ######### BBBBBBBBBBBBBBBBBBBBBBBBBB #########
    #########
    ######### STEP updating info with new data point y_{t+1}
    #########

    #########
    #########
    ######### 1 ######### updating 3-point already in info
    #########
    #########
    for(l in 1:nrows[t+1]) # update m with new index (time step) t+1
    {
      k <- info$k[l]
      j <- info$j[l]
      i <- info$i[l]

      if(is.na(j))  #update surface later
      {
        res <- eval2D_q_0(costQ, cumy1, cumy2, cumyS, k, t+1, beta)
        info$p1[l] <- res$p[1]
        info$p2[l] <- res$p[2]
        info$m[l] <- res$m
      }
      if((!is.na(j)) & is.na(i)) #update surface later
      {
        res <- eval2D_q_1(costQ, cumy1, cumy2, cumyS, j, k, t+1, beta)
        info$p1[l] <- res$p[1]
        info$p2[l] <- res$p[2]
        info$m[l] <- res$m
      }
      if((!is.na(j)) & (!is.na(i))) ### update only the argmin
      {
        if(info$surface[l] == 1)
        {
          res <- eval2D_q_21(costQ, cumy1, cumy2, cumyS, i, j, k, t+1, beta)
          info$m[l] <- res$m
          ### ### ### ### ### testing against new t+1 (delete later) ### ### ### ### ###
          testValue <- eval2D_q(costQ, cumy1, cumy2, cumyS, t+1, t+1, beta, info$p1[l], info$p2[l])
          if(testValue < info$m[l]){info$surface[l] <- 0} ### to remove in cleaning step
        }
        if(info$surface[l] == 2)
        {
          res <- eval2D_q_22(costQ, cumy1, cumy2, cumyS, i, j, k, t+1, beta)
          info$m[l] <- res$m
          ### ### ### ### ### testing against new t+1 (delete later) ### ### ### ### ###
          testValue <- eval2D_q(costQ, cumy1, cumy2, cumyS, t+1, t+1, beta, info$p1[l], info$p2[l])
          if(testValue < info$m[l]){info$surface[l] <- 0} ###  to remove in cleaning step
        }
      }
    }

    #########
    #########
    ######### 2 ######### adding new 3-points with t+1
    #########
    #########

    ### ### ### new quadratics in t+1 (always on the surface)

    temp <- eval2D_q_0(costQ, cumy1, cumy2, cumyS, t+1, t+1, beta)
    info <- rbind(info, c(t+1, NA, NA, temp$m, temp$p[1], temp$p[2], 1)) #update surface later

    ### ### ### new triple intersection points with  t+1 (always on the surface)
    ### ### ### new triple intersection points with  t+1 (always on the surface) CAN BE DONE BETTER...
    ### ### ### new triple intersection points with  t+1 (always on the surface)
    seen_indices <- NULL
    for(j in indexSet)
    {
      for(i in indexSet)
      {
        if((i < j) & eval2D_circleIntersection_ijk(costQ, cumy1, cumy2, cumyS, i, j, t+1))
        {
          temp1 <- eval2D_q_21(costQ, cumy1, cumy2, cumyS, i, j, t+1, t+1, beta)
          temp2 <- eval2D_q_22(costQ, cumy1, cumy2, cumyS, i, j, t+1, t+1, beta)
          surface1 <- 1
          surface2 <- 2
          ########## COVERING TEST ##########
          for(l in indexSet)
          {
            if((l != i) & (l != j))
            {
              testValue <- eval2D_q(costQ, cumy1, cumy2, cumyS, l, t+1, beta, temp1$p[1], temp1$p[2])
              if (testValue < temp1$m){surface1 <- 0}
              testValue <- eval2D_q(costQ, cumy1, cumy2, cumyS, l, t+1, beta, temp2$p[1], temp2$p[2])
              if (testValue < temp2$m){surface2 <- 0}
            }
          }
          if(surface1 == 1)
          {
            info <- rbind(info, c(t+1, j, i, temp1$m, temp1$p[1], temp1$p[2], 1))
            seen_indices <- unique(c(seen_indices, i, j))
          }
          if(surface2 == 2)
          {
            info <- rbind(info, c(t+1, j, i, temp2$m, temp2$p[1], temp2$p[2], 2))
            seen_indices <- unique(c(seen_indices, i, j))
          }
          #both  surface visible
        }
      }
    }

    ### ### ### new intersection i with t+1 (seen in triple ij(t+1) points on the surface)
    for(i in seen_indices)
    {
      temp <- eval2D_q_1(costQ, cumy1, cumy2, cumyS, i, t+1, t+1, beta)
      info <- rbind(info, c(t+1, i, NA, temp$m, temp$p[1], temp$p[2], 1)) #update surface later
    }

    ### ### ### unseen indices really unseen ?
    un_seen_indices <- setdiff(indexSet, seen_indices)
    #for(i in un_seen_indices)
    for(i in indexSet)
    {
      isOnSurface <- TRUE
      temp <- eval2D_q_1(costQ, cumy1, cumy2, cumyS, i, t+1, t+1, beta)

      for(j in indexSet)
      {
        if(j != i)
        {
          testValue <- eval2D_q(costQ, cumy1, cumy2, cumyS, j, t+1, beta, temp$p[1], temp$p[2])
          if(testValue < temp$m) # i is hidden by j
          {
            isOnSurface <- FALSE
            break
          }
        }
      }
      if(isOnSurface == TRUE)
      {
        info <- rbind(info, c(t+1, i, NA, temp$m, temp$p[1], temp$p[2], 1)) #update surface later
      }
    }

    #####
    ##### cleaning
    #####
    #UnseentriplePoints <- (!is.na(info$j)) & (!is.na(info$i)) & (info$surface == 0)
    #info <- info[!UnseentriplePoints,]

    ### add new index t
    indexSet <- c(indexSet, t+1)


    #########
    #########
    ######### 3 ######### UPDATE status SURFACE !!! for all q_0 and q_1
    #########
    #########

    nrows[t+1] <- nrow(info) ### count number of rows in info data-frame

    for(l in 1:nrows[t+1])
    {
      k <- info$k[l]
      j <- info$j[l]
      i <- info$i[l]
      ### ### ### ### ### ### ### ### ### ###
      if((is.na(j)) & (is.na(i)))
      {
        seen_surface <- 1
        for(r in indexSet)
        {
          if(k != r)
          {
            testValue <- eval2D_q(costQ, cumy1, cumy2, cumyS, r, t+1, beta, info$p1[l], info$p2[l])
            if(testValue < info$m[l]){seen_surface <- 0; break}
          }
        }
        if(seen_surface == 0){info$surface[l] <- 0}else{info$surface[l] <- 1}
      }
      ### ### ### ### ### ### ### ### ### ###
      if((!is.na(j)) & (is.na(i)))
      {
        seen_surface <- 1
        for(r in indexSet)
        {
          if((k != r) & (j != r))
          {
            testValue <- eval2D_q(costQ, cumy1, cumy2, cumyS, r, t+1, beta, info$p1[l], info$p2[l])
            if(testValue < info$m[l]){seen_surface <- 0; break}
          }
        }
        if(seen_surface == 0){info$surface[l] <- 0}else{info$surface[l] <- 1}
      }
      ### ### ### ### ### ### ### ### ### ###
    }

    nb[t+1] <- length(indexSet) ### count number of elements for DP
    #rows[t+1] <- nrow(info) ### count number of rows in info data-frame

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


