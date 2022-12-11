
coeff <- c(5,2,6,3,-1,-10)
E1 <- points.on.ellipse(coeff)
res <- ellipseRotation(coeff,NULL)
newcoeff <- res$coeff
E2 <- points.on.ellipse(newcoeff)

res$coeff
res$angle
newcoeff2 <- ellipseRotation(newcoeff, res$angle)


xlim = c(min(E1[,1],E2[,1]),max(E1[,1],E2[,1]))
ylim = c(min(E1[,2],E2[,2]),max(E1[,2],E2[,2]))


E4 <- points.on.ellipse(newcoeff2$coeff, 1000)



base::plot(E1, col = 1, type = 'l', asp = 1,
           xlim = xlim,
           ylim = ylim)
par(new = TRUE)
base::plot(E2, type = 'l', asp = 1, col = "green",
           xlim = xlim,
           ylim = ylim)
par(new = TRUE)
base::plot(E4, type = 'l', asp = 1, col = "red",
           xlim = xlim,
           ylim = ylim)
points(0,0, col = 2)
points(ellipseCenter(coeff), col = 1)


##############################################################################


#par(new = TRUE)
#base::plot(E3, type = 'l', asp = 1, col = "red",
#           xlim = xlim,
#           ylim = ylim)


#####

coeff1 <- c(10,7,6,1,-1,-3)

rotation <- res$rot
#rotation[1,1] <- abs(rotation[1,1])
#rotation[2,2] <- rotation[2,2]
#rotation[1,2] <- -rotation[1,2]
res2 <- ellipseRotation(coeff1, rotation)

coeff2 <- res2$coeff

Enew1 <- my.ellipse(coeff1)
Enew2 <- my.ellipse(coeff2)

base::plot(Enew1, col = 4, type = 'l', asp = 1,
           xlim = xlim,
           ylim = ylim)
par(new = TRUE)
base::plot(Enew2, type = 'l', asp = 1, col = 5,
           xlim = xlim,
           ylim = ylim)
points(0,0, col = 2)









