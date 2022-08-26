

###########################################################
#############     theta_grid_square_2D    #################
###########################################################

theta_grid_square_2D <- function(data, bound1 = c(0,1), bound2 = c(0,1), n12 = c(100,100), margin = c(0.05, 0.05))
{
    start1 = min(data$y1) - margin[1]*(max(data$y1) - min(data$y1))
    end1 = max(data$y1) + margin[1]*(max(data$y1) - min(data$y1))
    start2 = min(data$y2) - margin[2]*(max(data$y2) - min(data$y2))
    end2 = max(data$y2) + margin[2]*(max(data$y2) - min(data$y2))

    len <- max(end1-start1, end2-start2)
    m <- c(mean(data$y1), mean(data$y2))
    theta1 <- seq(m[1] - len/2, m[2] + len/2, length.out = n12[1])
    theta2 <- rev(seq(m[1] - len/2, m[2] + len/2, length.out = n12[2]))

    w <- expand.grid(theta2, theta1)
    Theta <- as.matrix(w)
    return(list(theta1 = theta1, theta2 = theta2, Theta = Theta))
}


rotate <- function(x) t(apply(x, 2, rev))


##############################################
#############     my_plot    #################
##############################################


my_plot <- function(M, iter, z, palette, legend = FALSE)
{
  par(mar = c(3, 3, 5, 6), xpd = TRUE)
  image(rotate(M),
        col= palette,
        zlim=c(0, z),
        xaxt= "n", yaxt= "n",
        main=paste('Iteration: n = ', as.character(iter)),
        asp = 1)
  if(legend == TRUE)
  {
    legend("topright",
           inset = c(-0.17, 0.0),
           legend = c(1:z),
           fill = palette,
           border = "black")
  }
}

##############################################
#############     plot_2D    #################
##############################################

plot_2D <- function(data, theta, beta = 4*length(data$y1), color = "ramp")
{
  #########
  ###
  ### DATA preprocessing
  ###
  y1 <- data$y1
  y2 <- data$y2
  theta1 <- theta$theta1
  theta2 <- theta$theta2
  Theta <- theta$Theta
  n <- length(y1)

  cumy1 <- cumsum(c(0, y1))
  cumy2 <- cumsum(c(0, y2))
  cumyS <- cumsum(c(0, y1^2 + y2^2))

  if(color == "grey"){palette <- gray.colors(n, start = 1, end = 0)}
  if(color == "ramp"){palette <- colorRampPalette(c("red", "white", "blue"), space = "Lab")(n)}

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
  eval_q_min <- function(k, t) ###minimum of q_{t}^{k}, data y_{k} to y_{t}
  {
    if(k == t){return(costQ[shift(k-1)] + beta)} ###costQ[shift(k)-1] = m_{k-1}
    return((t - k + 1)*eval_var(k,t) + costQ[shift(k-1)] + beta)
  }
  #########
  eval_q <- function(k, t, t1, t2) ###value of q_{t-1}^{k}(t1,t2)
  {
    return((t - k + 1)*((t1 - eval_meany1(k,t))^2 +  (t2 - eval_meany2(k,t))^2) + eval_q_min(k,t))
  }

  costQ <- rep(0, n + 1)
  costQ[1] <- -beta

  #########
  ###
  ### update rule
  ###
  M <- matrix(0, length(theta2), length(theta1))

  for(t in 1:n)
  {
    ###FIND min quadratics
    min_temp <- Inf
    for(k in 1:t)
    {
      val1 <- eval_q_min(k, t)
      if(val1 < min_temp){min_temp <- val1}
    }
    costQ[shift(t)] <- min_temp

    ###VALUES in Matrix
    for(s in 1:nrow(Theta))
    {
      min_temp <- Inf
      for(i in 1:t)
      {
        t_it1 <- Theta[s,1]
        t_it2 <- Theta[s,2]
        val2 <- eval_q(i, t, t_it1, t_it2)
        if(val2 < min_temp)
        {
          min_temp <- val2
          index <- i
        }
      }
      M[s] <- index
    }
    my_plot(M, t, n, palette, legend = TRUE)
  }
}



######################################################
#############     plot_Regression    #################
######################################################

plot_Regression <- function(data, theta, beta = 4*length(data$y1), color = "ramp")
{

}



