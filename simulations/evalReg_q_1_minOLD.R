
evalReg_q_1_minOLD <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t, beta)
{
  #print("evalReg_q_1_min")
  #print(c(j,k,t))
  coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
  coeffk <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
  coeff1C <- coeffj - coeffk

  newCoeff <- ellipseRotation(coeffk)
  coeffk <- newCoeff$coeff
  coeffj <- ellipseRotation(coeffj, newCoeff$rot)
  coeff1C <- ellipseRotation(coeff1C, newCoeff$rot)
  coeffj<- coeffj$coeff
  coeff1C<- coeff1C$coeff

  center <- ellipseCenter(coeff1C)
  centerj <- ellipseCenter(coeffj)
  centerk <- ellipseCenter(coeffk)


  if(ellipseEval(coeff1C, center$x, center$y) > 0){return(Inf)}

  print("TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTEST")
  print(center)
  print(centerk)
  print(centerj)

  fk <- function(s, type)
  {
    branch <- ellipseBranch(coeff1C, s, ref0 = centerk$y, type = type)
    return(ellipseEval(coeffk, branch$x, branch$y))
  }
  fj <- function(s, type)
  {
    branch <- ellipseBranch(coeff1C, s, ref0 = centerj$y, type = type)
    return(ellipseEval(coeffj, branch$x, branch$y))
  }

  #################################################################


  if((j < k) & (k <t) & (t = 7))
  {
    print("TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTEST")
    print(center)
    print(centerk)
    print(centerj)
    coeffj2 <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
    coeffk2 <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
    coeff1C2 <- coeffj2 - coeffk2

    newCoeff2 <- ellipseRotation(coeffk2)
    coeffk2 <- newCoeff2$coeff
    coeffj2 <- ellipseRotation(coeffj2, newCoeff2$rot)
    coeff1C2 <- ellipseRotation(coeff1C2, newCoeff2$rot)
    coeffj2<- coeffj2$coeff
    coeff1C2<- coeff1C2$coeff


    center2 <- ellipseCenter(coeff1C2)
    centerj2 <- ellipseCenter(coeffj2)
    centerk2 <- ellipseCenter(coeffk2)
    fk2 <- function(s, type)
    {
      branch <- ellipseBranch(coeff1C2, s, ref0 = centerk2$y,  type = type)
      return(evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta, branch$x, branch$y))
    }
    fj2 <- function(s, type)
    {
      branch <- ellipseBranch(coeff1C2, s, ref0 = centerj2$y, type = type)
      return(evalReg_q(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta, branch$x, branch$y))
    }
    nom <- 10000
    h <- 1/nom
    pl1 <- sapply(seq(-1,-h, length.out = nom), function(x) fk2(x, "left"))
    pl2 <- sapply(seq(0,1-h, length.out = nom), function(x) fk2(x, "left"))
    pl3 <- sapply(seq(1,h, length.out = nom), function(x) fk2(x, "right"))
    pl4 <- sapply(seq(0,-1+h, length.out = nom), function(x) fk2(x, "right"))

    par(mfrow=c(1,1))

    ########################################

    col <- c(rep(1,nom),rep(2,nom),rep(3,nom),rep(4,nom))

    plot(sqrt(c(pl1,pl2,pl3,pl4)),
         type = "p",
         pch = 20,
         lwd = 2,
         col = col)

    ########################################

    l1k <- golden_Search(fk2, -1, 0, "left")
    l2k <- golden_Search(fk2, 0, 1, "left")
    r1k <- golden_Search(fk2, -1, 0, "right")
    r2k <- golden_Search(fk2, 0, 1, "right")

    l1j <- golden_Search(fj2, -1, 0, "left")
    l2j <- golden_Search(fj2, 0, 1, "left")
    r1j <- golden_Search(fj2, -1, 0, "right")
    r2j <- golden_Search(fj2, 0, 1, "right")

    M <- c(l1k$m, l2k$m, r1k$m, r2k$m, l1j$m, l2j$m, r1j$m, r2j$m)
    ind <- which.min(M)
    print("COUTOUR ANALYSIS")
    print(M)
    test <- c(l1k$s, l2k$s, r1k$s, r2k$s, l1j$s, l2j$s, r1j$s, r2j$s)
    print(test)
    test <- (test < (-1+10^(-13))) + 2*(abs(test) < 10^(-13)) +  4*(test > (1-10^(-13)))
    print(test)
    test <- test + 8*(M == Inf)
    print(test + 8*(M == Inf))
    print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    print(sum(test == 0) == 4)
    print("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")

    print(isAnEllipse(coeff1C2))

    E1C <- my.ellipse(coeff1C2, n.points = 4*nom)
    cent <- ellipseCenter(coeffk2)
    GG <- ellipseEval(coeffk2, cent$x, cent$y)

    coeffk1 <- coeffk2
    coeffk2 <- coeffk2
    coeffk3 <- coeffk2
    coeffk4 <- coeffk2

    coeffk1[6] <-  coeffk2[6] - GG - 0.5
    coeffk2[6] <-  coeffk2[6] - GG  - 2
    coeffk3[6] <-  coeffk2[6] - GG -4
    coeffk4[6] <-  coeffk2[6] - GG  -6
    Ek1 <- my.ellipse(coeffk1)
    Ek2 <- my.ellipse(coeffk2)
    Ek3 <- my.ellipse(coeffk3)
    Ek4 <- my.ellipse(coeffk4)
    minx <- min(c(E1C[,1]))
    maxx <- max(c(E1C[,1]))
    miny <- min(c(E1C[,2]))
    maxy <- max(c(E1C[,2]))

    centerk2 <- ellipseCenter(coeffk2)

    print("test 2 test 2 test 2 test 2 test 2 test 2 ")
    print(centerk2)
    minEy <- which.min(E1C[,2])
    maxEy <- which.max(E1C[,2])

    xBAS <- E1C[minEy,1]
    xHAUT <- E1C[maxEy,1]
    yBAS <- min(E1C[,2])
    yHAUT <- max(E1C[,2])

    x <- E1C[,1]
    y <- E1C[,2]

    ligne <- (xHAUT - xBAS)*y-((yHAUT - yBAS)*x + yBAS*(xHAUT - xBAS)-(yHAUT - yBAS)*xBAS)

    Select1 <- (E1C[,2] <= centerk2$y) & (ligne >= 0)
    Select2 <- (E1C[,2] > centerk2$y) & (ligne >= 0)
    Select3 <- (E1C[,2] > centerk2$y) & (ligne < 0)
    Select4 <- (E1C[,2] <= centerk2$y) & (ligne < 0)
    E1C <- rbind(E1C[Select1,], E1C[Select2,], E1C[Select3,], E1C[Select4,])

    col <- c(rep(1,sum(Select1)),
             rep(2,sum(Select2)),
             rep(3,sum(Select3)),
             rep(4,sum(Select4)))


    base::plot(E1C, xlim = c(minx,maxx), ylim = c(miny,maxy),
               col = col, type = 'p', asp = 1)
    par(new = TRUE)
    base::plot(Ek1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek3, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek4, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 1, lwd = 3)
    points(centerk2$x, centerk2$y, pch = 16, col = colors()[459])
    print(coeffk1)
    A <- coeffk1[1]
    B <- coeffk1[2]
    C <- coeffk1[3]
    M <- matrix(c(A,B,B,C),2,2)
    sdvM <- svd(M)
    print(sdvM)
    segments(
      x0 = centerk2$x,
      y0 = centerk2$y,
      x1 = centerk2$x + sdvM$u[1,1],
      y1 = centerk2$y + sdvM$u[2,1], pch = 16, col = colors()[459]
    )
    segments(
      x0 = centerk2$x,
      y0 = centerk2$y,
      x1 = centerk2$x + sdvM$u[1,2],
      y1 = centerk2$y + sdvM$u[2,2], pch = 16, col = colors()[459]
    )

    #newCoeff <- ellipseRotation(coeffk1)
    #print(newCoeff)
    #Ekrot <- my.ellipse(newCoeff)
    #par(new = TRUE)
    #base::plot(Ekrot, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[422], type = 'l', asp = 1, lty = 3)



    stop("zgehterhrth ryjyrtj ryj ")
  }

  #################################################################


  # 2 explorations s in [0, 1] left and right
  # 2 explorations s in [-1, 0] left and right
  l1k <- golden_Search(fk, -1, 0, "left")
  l2k <- golden_Search(fk, 0, 1, "left")
  r1k <- golden_Search(fk, -1, 0, "right")
  r2k <- golden_Search(fk, 0, 1, "right")

  l1j <- golden_Search(fj, -1, 0, "left")
  l2j <- golden_Search(fj, 0, 1, "left")
  r1j <- golden_Search(fj, -1, 0, "right")
  r2j <- golden_Search(fj, 0, 1, "right")

  M <- c(l1k$m, l2k$m, r1k$m, r2k$m, l1j$m, l2j$m, r1j$m, r2j$m)
  ind <- which.min(M)
  print("COUTOUR ANALYSIS")
  print(M)
  test <- c(l1k$s, l2k$s, r1k$s, r2k$s, l1j$s, l2j$s, r1j$s, r2j$s)
  print(test)
  test <- (test < (-1+10^(-13))) + 2*(abs(test) < 10^(-13)) +  4*(test > (1-10^(-13)))
  print(test)
  test <- test + 8*(M == Inf)
  print(test + 8*(M == Inf))
  print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
  print(sum(test == 0) == 4)
  print("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")
  #print(ind)
  m <- M[ind]
  #s <- c(l1k$s, l2k$s, r1k$s, r2k$s, l1j$s, l2j$s, r1j$s, r2j$s)[ind]
  #type <- c("left", "left", "right", "right", "left", "left", "right", "right")[ind]
  return(m)
}





##############################################
#####
##### Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + G #####
#####
##### s in [-1,1] with y = ref0 obtained at s = 0
#####
##############################################

ellipseBranch <- function(coeff, s, ref0, type)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  ###
  ### y in ]yMin, yMax[ to have two roots in x (for ellipses)
  ###
  u <- A*C - B^2
  v <- B*D - A*E
  w <- v^2 + u*(D^2 - A*G)

  if((u < 0) | (w < 0)){return(NULL)}
  if((s < -1) | (s > 1)){return(NULL)}

  yMin <- (v-sqrt(w))/u
  yMax <- (v+sqrt(w))/u

  ###
  ### computing y value
  ###
  if(s >= 0){y <- ref0 + s*(yMax - ref0)}
  if(s < 0){y <- ref0 + s*(ref0 - yMin)}
  if(abs(s - 1)<10^(-14)){y <- yMax}
  if(abs(s + 1)<10^(-14)){y <- yMin}

  ###
  ### computing x value
  ###
  z1 <- B*y+D
  z2 <- z1^2 - A*(C*y^2 + 2*E*y + G)
  if((z2 < 10^(-13)) | (abs(s - 1)<10^(-14)) | (abs(s + 1)<10^(-14))){x <- -z1/A}
  else
  {
    if(type == "left"){x <- (-z1 - sqrt(z2))/A}
    if(type == "right"){x <- (-z1 + sqrt(z2))/A}
  }

  return(list(x = x, y = y))
}
