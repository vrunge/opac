
#############################################

## use after rotating coeff_eval
ellipseCoordinatesFromPoint <- function(coeff_constraint, point, theta, type)
{
  A <- coeff_constraint[1]
  B <- coeff_constraint[2]
  C <- coeff_constraint[3]
  D <- coeff_constraint[4]
  E <- coeff_constraint[5]
  G <- coeff_constraint[6]
  c1 <- point$x
  c2 <- point$y

  a <- A*cos(theta)^2 + 2*B*cos(theta)*sin(theta) + C*sin(theta)^2
  b <- D*cos(theta) + E*sin(theta) + A*c1*cos(theta) + B*c1*sin(theta) + B*c2*cos(theta) + C*c2*sin(theta)
  c <- A*c1^2 + 2*B*c1*c2 + C*c2^2 + 2*D*c1 + 2*E*c2 + G

  d <- b^2 - a*c
  if(d < 10^(-14)){e <- 0}else{e <- sqrt(d)/a}

  # center + (u*cos(theta), u*sin(theta))
  if(type == "plus")
  {
    x <- c1 + (-b/a + e)*cos(theta)
    y <- c2 + (-b/a + e)*sin(theta)
  }
  if(type == "minus")
  {
    x <- c1 + (-b/a - e)*cos(theta)
    y <- c2 + (-b/a - e)*sin(theta)
  }
  #print("test")
  #print(theta/(pi))
  #print(c(c1,c2))
  #print(c(-b/a - e))
  #print(c(-b/a + e))
  #print(e)
  #print(c(type))

  #if(theta/(pi) > 0.25)
  #{
  #  E <- points.on.ellipse(coeff_constraint, n.points = 1000) # POINTS ON THE ELLIPSE
  #  base::plot(E, col = colors()[459], type = 'l', asp = 1, lty = 1, lwd = 3)
  #  points(c1, c2, pch = 16, col = colors()[459])
  #  stop("stop")
  #}
  return(list(x = x, y = y))
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





##############################################
#####
##### Ax^2 + 2Bxy + Cy^2 + 2Dx + 2Ey + G #####
#####
##### s in [-1,1]. : type = left or right
#####
##############################################

ellipseBranch3 <- function(coeff, s, type)
{
  A <- coeff[1]
  B <- coeff[2]
  C <- coeff[3]
  D <- coeff[4]
  E <- coeff[5]
  G <- coeff[6]

  ##############################
  ###
  ### y in ]yMin, yMax[ to have two roots in x (for ellipses)
  ###
  uY <- A*C - B^2
  vY <- B*D - A*E
  wY <- vY^2 + uY*(D^2 - A*G)

  if((uY < 0) | (wY < 0)){return(NULL)}
  if((s < -1) | (s > 1)){return(NULL)}

  yMin <- (vY-sqrt(wY))/uY
  yMax <- (vY+sqrt(wY))/uY


  ##############################
  ###
  ### x in ]xMin, xMax[ to have two roots in y (for ellipses)
  ###
  uX <- A*C - B^2
  vX <- B*E - C*D
  wX <- vX^2 + uX*(E^2 - C*G)

  if((uX < 0) | (wX < 0)){return(NULL)}
  if((s < -1) | (s > 1)){return(NULL)}

  xMin <- (vX-sqrt(wX))/uX
  xMax <- (vX+sqrt(wX))/uX
  yLeft <- -(B*xMin + E)/C
  yRight <- -(B*xMax + E)/C

  ###
  ### computing y ref and s scaling
  ###
  if((type == "left")){ref0 <- yLeft}
  if((type == "right")){ref0 <- yRight}

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


## use after rotating by rotation of coeff_eval
evalReg_q_onConstraint2 <- function(coeff_eval, coeff_constraint, s, type)
{
  # find center
  point <- ellipseCenter(coeff_eval)
  # get point on quadrant (type = plus or minus, theta in [0, pi])
  pointToEval <- ellipseBranch(coeff_constraint, s, point$y, type)

  return(ellipseEval(coeff_eval, pointToEval$x, pointToEval$y))
}




## use after rotating by rotation of coeff_eval
evalReg_q_onConstraint3 <- function(coeff_eval, coeff_constraint, s, type)
{
  # get point on quadrant (type = plus or minus, s in [0, 1])
  pointToEval <- ellipseBranch3(coeff_constraint, s, type)

  return(ellipseEval(coeff_eval, pointToEval$x, pointToEval$y))
}






evalReg_q_1_minOLD2 <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t, beta)
{

  ### ### GET COEFF ### ###
  ### ### GET COEFF ### ###
  coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
  coeffk <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
  coeff1C <- coeffj - coeffk

  ########################################################################
  #if((j < k) & (k <t) & isAnEllipse(coeff1C)& (t == 7))
  #{
  # E3 <- points.on.ellipse(coeff1C, 1000)
  #  coeffk1 <- coeffk
  #  A <- coeffk[1]
  #  B <- coeffk[2]
  #  C <- coeffk[3]
  #  D <- coeffk[4]
  #  E <- coeffk[5]
  #  G <- coeffk[6]
  #  minim <-  (C*D^2 + A*E^2 - 2*B*E*D)/(A*C-B^2)
  #  coeffk1[6] <- minim - 6
  #  Ek1 <- points.on.ellipse(coeffk1, 1000)
  #  xlim = c(min(E3[,1]),max(E3[,1]))
  #  ylim = c(min(E3[,2]),max(E3[,2]))

  #  base::plot(E3, col = 1, type = 'l', asp = 1,
  #             xlim = xlim,
  #             ylim = ylim)
  #  par(new = TRUE)
  # base::plot(Ek1,  col = colors()[459], type = 'l', asp = 1, lty = 1, lwd = 3,
  #             xlim = xlim,
  #             ylim = ylim)
  #}
  ########################################################################

  ### ### ROTATION MATRIX ### ###
  ### ### ROTATION MATRIX ### ###
  temp <- ellipseRotation(coeffk)
  angle <- temp$angle #rotation MATRIX

  ### ### ROTATION ### ###
  ### ### ROTATION ### ###
  coeffk <- temp$coeff
  temp <- ellipseRotation(coeff1C, angle)
  coeff1C <- temp$coeff

  #temp <- ellipseRotation(coeffj, rot)
  #1coeffj <- temp$coeff


  ########################################################################
  if((j < k) & (k <t) & isAnEllipse(coeff1C) & (t == 7))
  {
    centerk <- ellipseCenter(coeffk)
    angle1 <- ellipseAngleFromConstant(coeff1C, centerk$x, "x")
    angle2 <- ellipseAngleFromConstant(coeff1C, centerk$y, "y")
    angle3 <- c(0, pi/2, pi, 3*pi/2)
    angles <- sort(c(angle1, angle2, angle3))

    if(ellipseEval(coeff1C, centerk$x, center$y) < 0){stop("center not in the constraint")}
    nom <- 10000
    h <- 1/nom
    pl <- sapply(seq(0, 2*pi - h, length.out = nom), function(x) evalReg_q_onConstraint(coeffk, coeff1C, x))

    par(mfrow=c(1,1))

    ################ plot value of k on the constraint (SQRT q_t^k) ###############

    col <- c(rep(1,nom), rep(2,nom), rep(3,nom), rep(4,nom))
    z <- sqrt(c(pl1,pl2,pl3,pl4))
    x <- seq(0,2*pi, length.out = length(z))

    plot(x[1:nom], z[1:nom], type = "l", main = "sqrt(value)", lwd = 7, col = 1,
         xlim = c(0,2*pi), ylim = c(min(z),max(z)))
    lines(x[(nom+1):(2*nom)],z[(nom+1):(2*nom)], col = 2, lwd = 7)
    lines(x[(2*nom+1):(3*nom)],z[(2*nom+1):(3*nom)], col = 3, lwd = 7)
    lines(x[(3*nom+1):(4*nom)],z[(3*nom+1):(4*nom)], col = 4, lwd = 7)

    ##############PLOT THE ELLIPSE CONSTRAINT #################
    centreK <- ellipseCenter(coeffk)
    coeffk1 <- coeffk; coeffk2 <- coeffk; coeffk3 <- coeffk; coeffk4 <- coeffk
    A <- coeffk[1]
    B <- coeffk[2]
    C <- coeffk[3]
    D <- coeffk[4]
    E <- coeffk[5]
    G <- coeffk[6]
    minim <-  (C*D^2 + A*E^2 - 2*B*E*D)/(A*C-B^2)
    coeffk1[6] <-  minim - 3
    coeffk2[6] <-  minim -6
    coeffk3[6] <-  minim -8
    coeffk4[6] <- minim - 10
    Ek1 <- my.ellipse(coeffk1)
    Ek2 <- my.ellipse(coeffk2)
    Ek3 <- my.ellipse(coeffk3)
    Ek4 <- my.ellipse(coeffk4)

    E1C <- my.ellipse(coeff1C, n.points = 4*nom) # POINTS ON THE ELLIPSE
    minx <- min(c(E1C[,1]))
    maxx <- max(c(E1C[,1]))
    miny <- min(c(E1C[,2]))
    maxy <- max(c(E1C[,2]))

    Select1 <- (E1C[,1] >= centreK$x) & (E1C[,2] >= centreK$y)
    Select2 <- (E1C[,1] < centreK$x) & (E1C[,2] >= centreK$y)
    Select3 <- (E1C[,1] < centreK$x) & (E1C[,2] < centreK$y)
    Select4 <- (E1C[,1] >= centreK$x) & (E1C[,2] < centreK$y)
    E1C <- rbind(E1C[Select1,], E1C[Select2,], E1C[Select3,], E1C[Select4,])
    col <- c(rep(1,sum(Select1)),
             rep(2,sum(Select2)),
             rep(3,sum(Select3)),
             rep(4,sum(Select4)))



    ####
    minEy <- which.min(E1C[,2])
    maxEy <- which.max(E1C[,2])

    xBAS <- E1C[minEy,1]
    xHAUT <- E1C[maxEy,1]
    yBAS <- min(E1C[,2])
    yHAUT <- max(E1C[,2])

    minEx <- which.min(E1C[,1])
    maxEx <- which.max(E1C[,1])

    yGAUCHE <- E1C[minEx,2]
    yDROITE <- E1C[maxEx,2]
    xGAUCHE <- min(E1C[,1])
    xDROITE <- max(E1C[,1])


    x <- E1C[,1]
    y <- E1C[,2]

    ligneY <- (xHAUT - xBAS)*y-((yHAUT - yBAS)*x + yBAS*(xHAUT - xBAS)-(yHAUT - yBAS)*xBAS)
    ligneX <- (xDROITE - xGAUCHE)*y-((yDROITE - yGAUCHE)*x + yGAUCHE*(xDROITE - xGAUCHE)-(yDROITE - yGAUCHE)*xGAUCHE)



    Select1 <- (ligneX <= 0) & (ligneY <= 0)
    Select2 <- (ligneX > 0) & (ligneY <= 0)
    Select3 <-  (ligneX > 0) & (ligneY > 0)
    Select4 <- (ligneX <= 0) & (ligneY > 0)
    E1C <- rbind(E1C[Select1,], E1C[Select2,], E1C[Select3,], E1C[Select4,])

    col <- c(rep(1,sum(Select1)),
             rep(2,sum(Select2)),
             rep(3,sum(Select3)),
             rep(4,sum(Select4)))
    ####



    base::plot(E1C, xlim = c(minx,maxx), ylim = c(miny,maxy),
               col = col, type = 'p', asp = 1)
    par(new = TRUE)
    base::plot(Ek1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek3, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek4, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    points(centreK$x, centreK$y, pch = 16, col = colors()[459])

    stop("STOP MON STOP")
  }
  ########################################################################


  # 2 explorations s in [0, 1] left and right
  # 2 explorations s in [-1, 0] left and right
  #l1k <- golden_Search(fk, -1, 0, "left")
  #l2k <- golden_Search(fk, 0, 1, "left")
  #r1k <- golden_Search(fk, -1, 0, "right")
  #r2k <- golden_Search(fk, 0, 1, "right")

  #l1j <- golden_Search(fj, -1, 0, "left")
  #l2j <- golden_Search(fj, 0, 1, "left")
  #r1j <- golden_Search(fj, -1, 0, "right")
  #r2j <- golden_Search(fj, 0, 1, "right")

  #M <- c(l1k$m, l2k$m, r1k$m, r2k$m, l1j$m, l2j$m, r1j$m, r2j$m)
  #ind <- which.min(M)
  #m <- M[ind]

  return(Inf)
}





evalReg_q_1_minOLD3 <- function(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, k, t, beta)
{

  ### ### GET COEFF ### ###
  ### ### GET COEFF ### ###
  coeffj <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, j, t, beta)
  coeffk <- ellipseCoeff(costQ, cumX, cumY, cumXY, cumSX, cumSY, k, t, beta)
  coeff1C <- coeffj - coeffk

  ########################################################################
  #if((j < k) & (k <t) & isAnEllipse(coeff1C)& (t == 7))
  #{
  # E3 <- points.on.ellipse(coeff1C, 1000)
  #  coeffk1 <- coeffk
  #  A <- coeffk[1]
  #  B <- coeffk[2]
  #  C <- coeffk[3]
  #  D <- coeffk[4]
  #  E <- coeffk[5]
  #  G <- coeffk[6]
  #  minim <-  (C*D^2 + A*E^2 - 2*B*E*D)/(A*C-B^2)
  #  coeffk1[6] <- minim - 6
  #  Ek1 <- points.on.ellipse(coeffk1, 1000)
  #  xlim = c(min(E3[,1]),max(E3[,1]))
  #  ylim = c(min(E3[,2]),max(E3[,2]))

  #  base::plot(E3, col = 1, type = 'l', asp = 1,
  #             xlim = xlim,
  #             ylim = ylim)
  #  par(new = TRUE)
  # base::plot(Ek1,  col = colors()[459], type = 'l', asp = 1, lty = 1, lwd = 3,
  #             xlim = xlim,
  #             ylim = ylim)
  #}
  ########################################################################

  ### ### ROTATION MATRIX ### ###
  ### ### ROTATION MATRIX ### ###
  temp <- ellipseRotation(coeffk)
  angle <- temp$angle #rotation MATRIX

  ### ### ROTATION ### ###
  ### ### ROTATION ### ###
  coeffk <- temp$coeff
  temp <- ellipseRotation(coeff1C, angle)
  coeff1C <- temp$coeff

  #temp <- ellipseRotation(coeffj, rot)
  #1coeffj <- temp$coeff


  ########################################################################
  if((j < k) & (k <t) & isAnEllipse(coeff1C) & (t == 7))
  {
    centerk <- ellipseCenter(coeffk)
    angle1 <- ellipseAngleFromConstant(coeff1C, centerk$x, "x")
    angle2 <- ellipseAngleFromConstant(coeff1C, centerk$y, "y")
    angle3 <- c(0, pi/2, pi, 3*pi/2)
    print(unlist(angle1)*180/pi)
    print(unlist(angle2)*180/pi)
    print(angle3*180/pi)
    angles <- sort(c(unlist(angle1), unlist(angle2), angle3))

    print(angles)
    #stop("ee")

    if(ellipseEval(coeff1C, centerk$x, center$y) < 0){stop("center not in the constraint")}
    nom <- 10000
    h <- 1/nom
    allAngles <- seq(0, 2*pi - h, length.out = nom)
    pl <- sapply(allAngles, function(x) evalReg_q_onConstraint(coeffk, coeff1C, x))

    par(mfrow=c(1,1))

    ################ plot value of k on the constraint (SQRT q_t^k) ###############
    l1 <- ((allAngles >= angles[1]) & (allAngles < angles[2]))
    l2 <- ((allAngles >= angles[2]) & (allAngles < angles[3]))
    l3 <- ((allAngles >= angles[3]) & (allAngles < angles[4]))
    l4 <- ((allAngles >= angles[4]) & (allAngles < angles[5]))
    l5 <- ((allAngles >= angles[5]) & (allAngles < angles[6]))
    l6 <- ((allAngles >= angles[6]) & (allAngles < angles[7]))
    l7 <- ((allAngles >= angles[7]) & (allAngles < angles[8]))
    l8 <- ((allAngles >= angles[8]))

    plot(allAngles[l1], pl[l1], type = "l", main = "sqrt(value)", lwd = 7, col = 1, xlim = c(0,2*pi), ylim = c(min(pl),max(pl)))
    lines(allAngles[l2], pl[l2], col = 2, lwd = 7)
    lines(allAngles[l3], pl[l3], col = 3, lwd = 7)
    lines(allAngles[l4], pl[l4], col = 4, lwd = 7)
    lines(allAngles[l5], pl[l5], col = 5, lwd = 7)
    lines(allAngles[l6], pl[l6], col = 6, lwd = 7)
    lines(allAngles[l7], pl[l7], col = 7, lwd = 7)
    lines(allAngles[l8], pl[l8], col = 8, lwd = 7)

    ##############PLOT THE ELLIPSE CONSTRAINT #################
    centreK <- ellipseCenter(coeffk)
    coeffk1 <- coeffk; coeffk2 <- coeffk; coeffk3 <- coeffk; coeffk4 <- coeffk
    A <- coeffk[1]
    B <- coeffk[2]
    C <- coeffk[3]
    D <- coeffk[4]
    E <- coeffk[5]
    G <- coeffk[6]
    minim <-  (C*D^2 + A*E^2 - 2*B*E*D)/(A*C-B^2)
    coeffk1[6] <-  minim - 3
    coeffk2[6] <-  minim -6
    coeffk3[6] <-  minim -8
    coeffk4[6] <- minim - 10
    Ek1 <- my.ellipse(coeffk1)
    Ek2 <- my.ellipse(coeffk2)
    Ek3 <- my.ellipse(coeffk3)
    Ek4 <- my.ellipse(coeffk4)

    E1C <- matrix(0,nom,2)

    for(i in 1:nom)
    {
      E1C[i,] <-  unlist(ellipseCoordinatesFromAngle(coeff1C, allAngles[i]))
    }

    minx <- min(c(E1C[,1]))
    maxx <- max(c(E1C[,1]))
    miny <- min(c(E1C[,2]))
    maxy <- max(c(E1C[,2]))
    print(angles)
    col <- (allAngles >= angles[1]) +
      (allAngles >= angles[2]) +
      (allAngles >= angles[3]) +
      (allAngles >= angles[4]) +
      (allAngles >= angles[5]) +
      (allAngles >= angles[6]) +
      (allAngles >= angles[7]) +
      (allAngles >= angles[8])

    base::plot(E1C, xlim = c(minx,maxx), ylim = c(miny,maxy), col = col, type = 'p', asp = 1)
    par(new = TRUE)
    base::plot(Ek1, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek2, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek3, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    par(new = TRUE)
    base::plot(Ek4, xlim = c(minx,maxx), ylim = c(miny,maxy), col = colors()[459], type = 'l', asp = 1, lty = 3)
    points(centreK$x, centreK$y, pch = 16, col = colors()[459])

    stop("STOP MON STOP")
  }
  ########################################################################


  # 2 explorations s in [0, 1] left and right
  # 2 explorations s in [-1, 0] left and right
  #l1k <- golden_Search(fk, -1, 0, "left")
  #l2k <- golden_Search(fk, 0, 1, "left")
  #r1k <- golden_Search(fk, -1, 0, "right")
  #r2k <- golden_Search(fk, 0, 1, "right")

  #l1j <- golden_Search(fj, -1, 0, "left")
  #l2j <- golden_Search(fj, 0, 1, "left")
  #r1j <- golden_Search(fj, -1, 0, "right")
  #r2j <- golden_Search(fj, 0, 1, "right")

  #M <- c(l1k$m, l2k$m, r1k$m, r2k$m, l1j$m, l2j$m, r1j$m, r2j$m)
  #ind <- which.min(M)
  #m <- M[ind]

  return(Inf)
}





