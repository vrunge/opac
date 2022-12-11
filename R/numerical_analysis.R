
########################################################################
########################################################################

golden_Search <- function(f, a, b)
{
  phi <- 2/(sqrt(5) + 1)

  x1 <- b - phi*(b - a)
  x2 <- a + phi*(b - a)

  f1 <- f(x1)
  f2 <- f(x2)

  if(f1 != Inf)
  {
    while (abs(f2 - f1) > 10^(-13))
    {
      if (f2 > f1)
      {
        b <- x2
        x2 <- x1
        f2 <- f1
        x1 <- b - phi*(b - a)
        f1 <- f(x1)
      }
      else
      {
        a <- x1
        x1 <- x2
        f1 <- f2
        x2 <- a + phi*(b - a)
        f2 <- f(x2)
      }
    }
  }
  ind <- which.min(c(f1, f2))

  return(list(s = c(x1, x2)[ind], m = c(f1, f2)[ind]))
}


########################################################################
########################################################################


Newton_Raphson <- function(f, fPrime, a, b, type)
{
  temp <- f(a)
  temp2 <- Inf
  if(type == "left"){c <- a}
  if(type == "right"){c <- b}

  while((abs(temp - temp2) > 10^(-13)) && (a <= c) && (c <= b))
  {
    temp <- f(c)
    c <- c - temp/fPrime(c)
    temp2 <- f(c)
  }
  if(c < a){return(NULL)}
  if(b < c){return(NULL)}
  return(c)
}








