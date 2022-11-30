

data <- dataGenerator_2D(chpts = c(30,100,120), means1 = c(0,5,0), means2 = c(7,1,-4))

4*log(120)

OP_2D(data, beta = 1)
OP_2D(data, beta = 7)
OP_2D(data, beta = 20)


OP_2D_PELT(data, beta = 1)
OP_2D_PELT(data, beta = 7)
OP_2D_PELT(data, beta = 100)



data <- dataGenerator_2D(chpts = 100, means1 = 0, means2 = 0)
res <- OP_2D_PELT(data, beta = 100)
res$nb
diff(res$nb)



####
n <- 2000
data <- dataGenerator_2D(n)

system.time(OP_2D_fast(data))
system.time(OP_2D(data))


OP_2D_fast(data)
OP_2D(data)



