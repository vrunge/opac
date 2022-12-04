


library(lamW)


#x > -1/e

y <- lambertW0(-exp(-1))
y

y*exp(y)


poisson <- function(A,B,C,t)
{
  return(A*t - B*log(t) + C)
}

A <- 10
B <-  1
C <- -5
A
B
-A/B * exp(C/B)
lambertW0(-A/B * exp(C/B))
lambertWm1(-A/B * exp(C/B))

result <- function(A,B,C)
{
  x <- -A/B * exp(C/B)
  if((x >= -exp(-1)) & x <= 0)return(c(t1 = (-B/A)*lambertW0(x),
                                          t2 =  (-B/A)*lambertWm1(x)))
  return(NULL)
}


poisson(A,B,C, result(A,B,C))



