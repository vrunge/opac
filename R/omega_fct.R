


omega_t_k_fct <- function(t, info, nrows,
                          indexSet, nb,
                          costQ, cumy1, cumy2, cumyS, beta)
{
  omega_t_k <- rep(-Inf, nb) #indexSet[u] = l, omega_t_k[u] for k = l

  for(i in 1:nrows)
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
      t1k <- eval_mean(cumy1, k, t)
      t2k <- eval_mean(cumy2, k, t)
      ind_k <- which(indexSet == k)
      # test min of q_t^jk invisible
      if(eval_q_min(costQ, cumy1, cumy2, cumyS, k, t, beta) > eval_q(costQ, cumy1, cumy2, cumyS, j, t, beta, t1k, t2k))
      {omega_t_k[ind_k] <- max(omega_t_k[ind_k], info$m[i])}
      t1j <- eval_mean(cumy1,j,t)
      t2j <- eval_mean(cumy2,j,t)
      ind_j <- which(indexSet == j)
      # test min of q_t^kj invisible
      if(eval_q_min(costQ, cumy1, cumy2, cumyS, j, t, beta) > eval_q(costQ, cumy1, cumy2, cumyS, k, t, beta, t1j, t2j))
      {omega_t_k[ind_j] <- max(omega_t_k[ind_j], info$m[i])}
    }
  }
  return(omega_t_k)
}

