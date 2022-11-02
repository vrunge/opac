




y1 <- data$y1
y2 <- data$y2
n <- length(y1)
cumy1 <- cumsum(c(0, y1))
cumy2 <- cumsum(c(0, y2))
cumyS <- cumsum(c(0, y1^2 + y2^2))

costQ <- rep(0, n + 1) # costQ[i] optimal cost for data y(1) to y(i-1)
costQ[1] <- -beta #costQ[2] = 0
