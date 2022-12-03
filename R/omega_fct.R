
###############################################
#############  omega_t_fct_1C #################
###############################################

omega_t_fct_1C <- function(t, info, nrows,
                             indexSet, nb,
                             costQ, cumy1, cumy2, cumyS, beta)
{
  omega_t <- rep(-Inf, nb)

  for(i in 1:nrows)
  {
    k <- info$k[i]
    j <- info$j[i]
    if(k == j)
    {
      ind <- which(indexSet == k)
      omega_t[ind] <- max(omega_t[ind], info$m[i])
    }
    else
    {
      #computing argmin of quadratics q_t^k (t1k,t2k)
      t1k <- eval_mean(cumy1, k, t)
      t2k <- eval_mean(cumy2, k, t)
      ind_k <- which(indexSet == k)
      # test min of q_t^jk invisible
      if(eval2D_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) > eval2D_q(costQ, cumy1, cumy2, cumyS, j, t, beta, t1k, t2k))
      {omega_t[ind_k] <- max(omega_t[ind_k], info$m[i])}

      #computing argmin of quadratics q_t^j (t1j,t2j)
      t1j <- eval_mean(cumy1, j, t)
      t2j <- eval_mean(cumy2, j, t)
      ind_j <- which(indexSet == j)
      #test min of q_t^kj invisible
      if(eval2D_q_min(costQ, cumy1, cumy2, cumyS, j, t, beta) > eval2D_q(costQ, cumy1, cumy2, cumyS, k, t, beta, t1j, t2j))
        {omega_t[ind_j] <- max(omega_t[ind_j], info$m[i])}
    }
  }
  return(omega_t)
}



###############################################
#############  omega_t_fct_2C #################
###############################################


omega_t_fct_2C <- function(t, info,
                          indexSet, nb,
                          costQ, cumy1, cumy2, cumyS, beta)
{
  omega_t <-  rep(Inf, nb)

  ####
  #### EXPLORE all integers in indexSet
  ####
  for(l in 1:nb) #omega_t[l]
  {
    #extract data
    k <- indexSet[l]
    selection <- apply(info[,1:3] == k, 1, function(x) any(x, na.rm = T))
    info_k <- info[selection , ]

    ###
    ### inside sub data-frame
    ###
    for(r in 1:nrow(info_k))
    {
      currentIndex <- c(info_k$i[r], info_k$j[r], info_k$k[r])
      #### m ####
      if(info_k$p1[r] != Inf)
      {
        p1 <- info_k$p1[r]
        p2 <- info_k$p2[r]
        m <- info_k$m[r]
        covered <- FALSE
        for(s in indexSet)
        {
          if(!(s %in% currentIndex))
          {
            temp <- eval2D_q(costQ, cumy1, cumy2, cumyS, s, t, beta, p1, p2)
            if(temp < m){covered <- TRUE; break}
          }
        }
        if(covered == FALSE){omega_t[l] <- min(omega_t[l], m)}
      }

      #### M ####
      if(info_k$q1[r] != Inf)
      {
        q1 <- info_k$q1[r]
        q2 <- info_k$q2[r]
        M <- info_k$M[r]
        covered <- FALSE
        for(s in indexSet)
        {
          if(!(s %in% currentIndex))
          {
            temp <- eval2D_q(costQ, cumy1, cumy2, cumyS, s, t, beta, q1, q2)
            if(temp < M){covered <- TRUE; break}
          }
        }
        if(covered == FALSE){omega_t[l] <- min(omega_t[l], M)}
      }
    }
  }
  return(omega_t)
}



################################################
#############   points_toInf   #################
################################################

points_toInf <- function(info, costQ, cumy1, cumy2, cumyS, t, beta)
{
  toTest <- which((info$k < t) & (!is.na(info$j)) & (!is.na(info$i)))

  for(l in toTest)
  {
    k <- info$k[l]
    j <- info$j[l]
    i <- info$i[l]
    ######
    if(info$p1[l] != Inf)
    {
      p1 <- info$p1[l]
      p2 <- info$p2[l]
      val_t_p <- eval2D_q(costQ, cumy1, cumy2, cumyS, t, t, beta, p1, p2)
      if(val_t_p < info$m[l])
      {
        info$p1[l] <- Inf
        info$p2[l] <- Inf
        info$m[l] <- Inf
      }
    }

    ######
    if(info$q1[l] != Inf)
    {
      q1 <- info$q1[l]
      q2 <- info$q2[l]
      val_t_p <- eval2D_q(costQ, cumy1, cumy2, cumyS, t, t, beta, q1, q2)
      if(val_t_p < info$M[l])
      {
        info$q1[l] <- Inf
        info$q2[l] <- Inf
        info$M[l] <- Inf
      }
    }
    ######
  }
  return(info)
}


