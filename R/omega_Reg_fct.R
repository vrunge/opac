
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
      # test min of q_t^jk invisible
      if(evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta) > evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1k, t2k))
      {omega_t[ind_k] <- max(omega_t[ind_k], info$m[i])}

      #computing argmin of quadratics q_t^j (t1j,t2j)
      coeff <- ellipseCenter(ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta))
      tj <- ellipseCenter(coeff)
      t1j <- tj$x
      t2j <- tj$y
      ind_j <- which(indexSet == j)
      #test min of q_t^kj invisible
      if(evalReg_q_min(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta) > evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, t1j, t2j))
      {omega_t[ind_j] <- max(omega_t[ind_j], info$m[i])}
    }
  }
  return(omega_t)
}

