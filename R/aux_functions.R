
#############################################
############     globalCost    ##############
#############################################

#' globalCost_2D
#' @description Computing the global cost of the segmentation for bivariate independent problem
#' @param data a dataframe with two components: y1 and y2, time series of same length
#' @param chpts vector of change-points obtained by one of the OP algorithm (2D case) using data-points 'data'
#' @param beta penalty value
#' @return the global cost value (minimal value for the optimization problem)
#' @examples
#' data <- dataGenerator2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4))
#' res <- OP_2D(data)
#' globalCost_2D(data, res$changepoints, 4*log(nrow(data)))
globalCost_2D <- function(data, chpts, beta)
{
  y1 <- data$y1
  y2 <- data$y2
  res <- 0
  chpts <- c(0, chpts)
  for(i in 1:(length(chpts)-1))
  {
    res <- res + sum((y1[(chpts[i]+1):chpts[i+1]] - mean(y1[(chpts[i]+1):chpts[i+1]]))^2 + (y2[(chpts[i]+1):chpts[i+1]] - mean(y2[(chpts[i]+1):chpts[i+1]]))^2)
  }
  res <- res + (length(chpts) - 2) * beta
  return(res)
}



#' globalCost_Reg
#' @description Computing the global cost of the segmentation for linear regression problem
#' @param data a dataframe with two components: x and y, time series of same length
#' @param chpts vector of change-points obtained by one of the OP algorithm (Reg case) using data-points 'data'
#' @param beta penalty value
#' @return the global cost value (minimal value for the optimization problem)
#' @examples
#' data <- dataGeneratorRegression(chpts = c(40,90), A = c(2,-1),  B = c(-1,2), meansX = c(1,2))
#' res <- OP_Reg(data)
#' globalCost_Reg(data, res$changepoints, 4*log(nrow(data)))
globalCost_Reg <- function(data, chpts, beta)
{
  x <- data$x
  y <- data$y
  res <- 0
  chpts <- c(0, chpts)
  for(i in 1:(length(chpts)-1))
  {
    xt <- x[(chpts[i]+1):chpts[i+1]]
    yt <- y[(chpts[i]+1):chpts[i+1]]
    res_lm <- lm(yt~xt)
    a <- unname(res_lm$coefficients[2])
    b <- unname(res_lm$coefficients[1])
    res <- res + sum((y[(chpts[i]+1):chpts[i+1]] - (a*x[(chpts[i]+1):chpts[i+1]] + b))^2)
  }
  res <- res + (length(chpts) - 2) * beta
  return(res)
}


