# Test functions

#---------------------------- burt -------------------------------------

data_burt <- function(n = 256, range = 1, rsnr = 4){
  
  sample_x <- seq(0, range, length.out = n)
  true_y <- 20*sample_x*cos(16*sample_x^(1.2)) - 20*as.numeric(identity(sample_x < 0.5))
  
  sdf <- sd(true_y)
  sigma <- sdf/rsnr
  
  epsilon <- rnorm(n, 0, sigma)
  
  sample_y <- true_y + epsilon
  
  position = which(as.numeric(identity(sample_x < 0.5) == 0) > 0)[1] - 1
  true_CP <- sample_x[position]
  
  data = list("true_y" = true_y, "sample_x" = sample_x, "sample_y" = sample_y, "true CPs" = true_CP, "noise" = sigma)
  return(data)
  
}

#---------------------------- cosine -------------------------------------

data_cosine <- function(n = 256, range = 1, rsnr = 4){
  
  sample_x <- seq(0, range, length.out = n)
  true_y <- cos(5.5*pi*sample_x) - 4*sign(0.23-sample_x) - 2*sign(0.3-sample_x) - 1.75*sign(0.55-sample_x) + 3*sign(0.7-sample_x)
  
  sdf <- sd(true_y)
  sigma <- sdf/rsnr
  
  epsilon <- rnorm(n, 0, sigma)
  
  sample_y <- true_y + epsilon
  
  position1 = which(sign(0.23-sample_x) == -1)[1] - 1
  position2 = which(sign(0.3-sample_x) == -1)[1] - 1
  position3 = which(sign(0.55-sample_x) == -1)[1] - 1
  position4 = which(sign(0.7-sample_x) == -1)[1] - 1
  position = c(position1, position2, position3, position4)
  
  true_CP <- sample_x[position]
  
  data = list("true_y" = true_y, "sample_x" = sample_x, "sample_y" = sample_y, "true CPs" = true_CP, "noise" = sigma)
  return(data)
  
}

#---------------------------- heavisine -------------------------------------

data_heavisine <- function(n = 256, range = 1, rsnr = 4){
  
  sample_x <- seq(0, range, length.out = n)
  true_y <- 4*sin(4*pi*sample_x) - sign(sample_x-0.3) - sign(0.72-sample_x)
  
  sdf <- sd(true_y)
  sigma <- sdf/rsnr
  
  epsilon <- rnorm(n, 0, sigma)
  
  sample_y <- true_y + epsilon
  
  position1 = which(sign(sample_x-0.3) == 1)[1] - 1
  position2 = which(sign(0.72-sample_x) == -1)[1] - 1
  position = c(position1, position2)
  
  true_CP <- sample_x[position]
  
  data = list("true_y" = true_y, "sample_x" = sample_x, "sample_y" = sample_y, "true CPs" = true_CP, "noise" = sigma)
  return(data)
  
  
}

#---------------------------- blip -------------------------------------

data_blip <- function(n = 256, range = 1, rsnr = 4){
  
  sample_x <- seq(0, range, length.out = n)
  true_y <- (0.32 + 0.6*sample_x + 0.3*exp(-100*(sample_x - 0.3)^2))*as.numeric(identity(sample_x <= 0.8)) + (-0.28 + 0.6*sample_x + 0.3*exp(-100*(sample_x - 1.3)^2))*as.numeric(identity(sample_x > 0.8)) 
  sdf <- sd(true_y)
  sigma <- sdf/rsnr
  
  epsilon <- rnorm(n, 0, sigma)
  
  sample_y <- true_y + epsilon
  
  position = which(as.numeric(identity(sample_x <= 0.8) == 0) > 0)[1] - 1
  
  true_CP <- sample_x[position]
  
  data = list("true_y" = true_y, "sample_x" = sample_x, "sample_y" = sample_y, "true CPs" = true_CP, "noise" = sigma)
  return(data)
  
  
}

#---------------------------- cubic -------------------------------------

data_poly3 <- function(n = 256, range = 6, rsnr = 4, gap1 = 5, gap2 = -8){ 
  
  t <- seq(-range/2, range/2, length.out = n)
  n1 <- as.integer(n/3)
  n2 <- n - 2*n1
  t1_sample <- t[1:n1]
  t2_sample <- t[(n1+1):(2*n1)]
  t3_sample <- t[-(1:(2*n1))]
  sample_x <- c(t1_sample, t2_sample, t3_sample)
  
  y1_true <- 2*t1_sample^3 + 3*t1_sample^2 - 0.6*t1_sample + gap1
  y2_true <- 2*t2_sample^3 + 3*t2_sample^2 - 0.6*t2_sample + gap2
  y3_true <- 2*t3_sample^3 + 3*t3_sample^2 - 0.6*t3_sample
  
  true_y <- c(y1_true, y2_true, y3_true)
  
  sdf <- sd(true_y)
  sigma <- sdf/rsnr
  
  epsilon <-  rnorm(n, 0, sigma)
  
  sample_y <- true_y + epsilon
  
  sample_x <- seq(0, range, length.out = n)
  
  true_CP = sample_x[c(n1, 2*n1)]
  data = list("true_y"=true_y, "sample_x"=sample_x, "sample_y"=sample_y, "true CPs"=true_CP)
  return(data)
  
}

#---------------------------- step -------------------------------------

data_pc2 <- function(n = 256, range = 1, rsnr = 4){ 
  
  t <- seq(0, 1, length.out = n)
  n1 <- as.integer(n/3)
  n2 <- n - 2*n1
  t1_sample <- t[1:n1]
  t2_sample <- t[(n1+1):(2*n1)]
  t3_sample <- t[-(1:(2*n1))]
  sample_x <- c(t1_sample, t2_sample, t3_sample)
  
  y1_true <- rep(3, length(t1_sample))
  y2_true <- rep(5, length(t2_sample))
  y3_true <- rep(1, length(t3_sample))
  
  true_y <- c(y1_true, y2_true, y3_true)
  
  sdf <- sd(true_y)
  sigma <- sdf/rsnr
  
  epsilon <-  rnorm(n, 0, sigma)
  
  sample_y <- true_y + epsilon
  
  sample_x <- seq(0, range, length.out = n)
  
  true_CP = sample_x[c(n1, 2*n1)]
  data = list("true_y"=true_y, "sample_x"=sample_x, "sample_y"=sample_y, "true CPs"=true_CP)
  return(data)
  
}

