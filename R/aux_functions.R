
#############################################
############     globalCost    ##############
#############################################

#' globalCost
#' @description Computing the global cost
#' @param data a dataframe with two components: y1 and y2, time series of same length
#'@param chpts vector of change-points obtained by one of the OP algorithl using data-points 'data'
#' @param beta penalty
#' @return the global cost value (minimal value for the optimization problem)
#' @examples
#' data <- dataGenerator2D(chpts = c(30,100,120), means1 = c(0,1,0), means2 = c(7,1,-4))
#' res <- OP_2D(data)
#' globalCost(data, res$changepoints, 4*log(nrow(data)))
globalCost <- function(data, chpts, beta)
{
  y1 <- data$y1
  y2 <- data$y2
  res <- 0
  chpts <- c(0, chpts)
  for(i in 1:(length(chpts)-1))
  {
    res <- res + sum((y1[(chpts[i]+1):chpts[i+1]] - mean(y1[(chpts[i]+1):chpts[i+1]]))^2 + (y2[(chpts[i]+1):chpts[i+1]] - mean(y2[(chpts[i]+1):chpts[i+1]]))^2)
  }
  res <- res + (length(chpts)-2)*beta
  return(res)
}
