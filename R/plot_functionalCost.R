
########################################################
#############     theta_grid_square    #################
########################################################


theta_grid_square <- function(data, bound1 = c(0,1), bound2 = c(0,1), n12 = c(100,100), margin = c(0.05, 0.05))
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



var_pop <- function(x){mean((x - mean(x))^2)}
my_var <- function(y1, y2, i, j){return(var_pop(y1[i:j]) + var_pop(y2[i:j]))}
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
  y1 <- data$y1
  y2 <- data$y2
  theta1 <- theta$theta1
  theta2 <- theta$theta2
  Theta <- theta$Theta
  n <- length(y1)

  vect_m <- rep(0, n + 1)
  vect_m[1] <- -beta

  if(color == "grey")
  {
    palette <- gray.colors(n, start = 1, end = 0)
  }
  if(color == "ramp")
  {
    palette <- colorRampPalette(c("red", "white", "blue"), space = "Lab")(n)
  }

  for(t in 1:n)
  {
    M <- matrix(0, length(theta2), length(theta1))
    min_temp <- t * my_var(y1, y2, 1, t) + vect_m[1] + beta
    for(k in 2:t)
    {
      val1 <- (t - k + 1) * my_var(y1, y2, k, t) + vect_m[k] + beta
      if(val1 < min_temp){min_temp <- val1}
    }
    vect_m[t+1] <- min_temp
    index <- 0

    for(s in 1:nrow(Theta))
    {
      min_temp <- Inf
      for(i in 1:t)
      {
        val2 <- (t - i + 1)*((Theta[s,1] - mean(y1[i:t]))^2 + (Theta[s,2] - mean(y2[i:t]))^2) + (t - i + 1) * my_var(y1, y2, i, t) + vect_m[i] + beta
        if(val2 < min_temp){
          min_temp <- val2
          index <- i
          M[s] <- index
        }
      }
    }
    my_plot(M, t, n, palette, legend = TRUE)
  }
}



######################################################
#############     plot_Regression    #################
######################################################

plot_Regression <- function(data, theta, beta = 4*length(data$y1), color = "ramp")
{
  y1 <- data$y1
  y2 <- data$y2
  theta1 <- theta$theta1
  theta2 <- theta$theta2
  Theta <- theta$Theta
  n <- length(y1)

  vect_m <- rep(0, n + 1)
  vect_m[1] <- -beta

  if(color == "grey")
  {
    palette <- gray.colors(n, start = 1, end = 0)
  }
  if(color == "ramp")
  {
    palette <- colorRampPalette(c("red", "white", "blue"), space = "Lab")(n)
  }

  for(t in 1:n)
  {
    M <- matrix(0, length(theta2), length(theta1))
    min_temp <- t * my_var(y1, y2, 1, t) + vect_m[1] + beta
    for(i in 2:t)
    {
      val1 <- (t - i + 1) * my_var(y1, y2, i, t) + vect_m[i] + beta
      if(val1 < min_temp){min_temp <- val1}
    }
    vect_m[t+1] <- min_temp
    index <- 0

    for(s in 1:nrow(Theta))
    {
      min_temp <- Inf
      for(i in 1:t)
      {
        val2 <- (t - i + 1)*(mean(y1[i:t]^2)*Theta[s,1]^2 + Theta[s,2]^2 + 2*mean(y1[i:t])*Theta[s,1]*Theta[s,2] - 2*mean(y2[i:t])*Theta[s,2] - 2*mean(y1[i:t]*y2[i:t])*Theta[s,1] + mean(y2[i:t]^2)) + vect_m[i] + beta
        if(val2 < min_temp){
          min_temp <- val2
          index <- i
          M[s] <- index
        }
      }
    }
    my_plot(M, t, n, palette, legend = TRUE)
  }
}



