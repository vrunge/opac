

# Plot Gauss data

n <- 15
mean <- 0
var <- 1
beta <- 4*log(n) # 4 better than 2

data <- dataGenerator2D(n)
plot(data$y1, data$y2)

theta <- theta_grid_square(data = data, n12 = c(100,100))
plot(theta$Theta)

plot_2D(data, theta, beta = beta, color = "ramp")
