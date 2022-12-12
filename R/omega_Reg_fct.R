
###############################################
#############  omega_t_fct_1C #################
###############################################

omega_t_Reg_fct_1C <- function(t, info, nrows,
                           indexSet, nb,
                           costQ, cumX, cumY, cumXY, cumSX, cumSY, beta)
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
      tk <- ellipseCenter(ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta))
      t1k <- tk$x
      t2k <- tk$y
      ind_k <- which(indexSet == k)
      # test min of q_t^k invisible
      if(evalReg_q_min2(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta) > evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta, t1k, t2k))
        {omega_t[ind_k] <- max(omega_t[ind_k], info$m[i])}

      #computing argmin of quadratics q_t^j (t1j,t2j)
      tj <- ellipseCenter(ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta))
      t1j <- tj$x
      t2j <- tj$y
      ind_j <- which(indexSet == j)
      #test min of q_t^j invisible
      if(evalReg_q_min2(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta) > evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1j, t2j))
        {omega_t[ind_j] <- max(omega_t[ind_j], info$m[i])}
    }
  }
  return(omega_t)
}



##########################################################
#############  omega_t_Reg_fct_1C_Approx #################
##########################################################

omega_t_Reg_fct_1C_Approx <- function(t, info, nrows,
                                      indexSet, nb,
                                      costQ, cumX, cumY, cumXY, cumSX, cumSY, beta)
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
      tk <- ellipseCenter(ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta))
      t1k <- tk$x
      t2k <- tk$y
      ind_k <- which(indexSet == k)
      # test min of q_t^k invisible
      #testApprox <- evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta, t1k, t2k)
      coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
      centerj <- ellipseCenter(coeffj)
      t1jt <- centerj$x
      t2jt <- centerj$y
      Dj2 <- (t1jt - t1k)^2 + (t2jt - t2k)^2
      A <- coeffj[1]
      B <- coeffj[2]
      C <- coeffj[3]
      Rj <- (1/2)*(A+C - sqrt((A-C)^2 + 4*B^2))
      testApprox <- Rj*Dj2 + evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
      #testApprox <- evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta, t1k, t2k)

      if(evalReg_q_min2(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta) > testApprox)
      {omega_t[ind_k] <- max(omega_t[ind_k], info$m[i])}

      #computing argmin of quadratics q_t^j (t1j,t2j)
      tj <- ellipseCenter(ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta))
      t1j <- tj$x
      t2j <- tj$y
      ind_j <- which(indexSet == j)
      #test min of q_t^j invisible
      #testApprox <- evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1j, t2j)
      coeffk <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
      centerk <- ellipseCenter(coeffk)
      t1kt <- centerk$x
      t2kt <- centerk$y
      Dk2 <- (t1kt - t1j)^2 + (t2kt - t2j)^2
      A <- coeffk[1]
      B <- coeffk[2]
      C <- coeffk[3]
      Rk <- (1/2)*(A+C - sqrt((A-C)^2 + 4*B^2))
      testApprox <- Rk*Dk2 +  evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
      #testApprox <- evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1j, t2j)

      if(evalReg_q_min2(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta) > testApprox)
      {omega_t[ind_j] <- max(omega_t[ind_j], info$m[i])}
    }
  }
  return(omega_t)
}

