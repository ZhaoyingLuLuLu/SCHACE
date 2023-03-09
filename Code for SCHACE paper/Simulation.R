library(devtools)
install_github("ZhaoyingLuLuLu/SCHACE/SCHACE-v2")
library(SCHACE)

### 3-fold cross-validation procedure with trimmed L2-error (percent = 10%)

## 100 replications

replications <- 100

CP <- rep(0, replications)
TP <- rep(0, replications)
MSE <- rep(0, replications)
FDR <- rep(0, replications)
df <- rep(0, replications)
lambda <- rep(0, replications)

for (i in 1:replications){

  set.seed(i)
  data_sample <- data_burt(rsnr = 4)
  sample_x <- data_sample$"sample_x"
  y <- data_sample$"sample_y"
  true_y <- data_sample$"true_y"
  true_CP <- data_sample$"true CPs"
  h <- 0.01
  range <- 1

  results_CV <- main.SCHACE(
    y,
    folds = 3,
    Bdf_set = seq(3, 15, 1),
    tuning = "crossvalidation",
    methods = "Trim",
    clambda = "lambdamin",
    percent = 0.1
  )

  df[i] <- results_CV$df
  lambda[i] <- results_CV$lambda

  CP[i] <- results_CV$`number of CPs`

  infor.result <- Informative(sample_x[results_CV$`locations of detected CPs`], true_CP, range)
  Informative_Sets <- infor.result$"infor"
  Uninformative_Sets <- infor.result$"uninfor"

  catch_results <- catch(sample_x, Informative_Sets, true_CP, h)
  TP[i] <- catch_results$"Number_of_Catched_CPs"

  FDR[i] <- (CP[i]-TP[i])/max(CP[i], 1)
  TDR <- TP[i] / length(true_CP)

  MSE[i] <-  mean((results_CV$`predicted y` - true_y)^2)
  
  print(i)

}


### extended Bayesian Information Criteria

## 100 replications

replications <- 100

df <- rep(0, replications)
lambda <- rep(0, replications)
CP <- rep(0, 100)
TP <- rep(0, 100)
MSE <- rep(0, 100)
FDR <- rep(0, 100)

for (i in 1:replications){

  set.seed(i)
  data_sample <- data_burt(rsnr = 4)
  sample_x <- data_sample$"sample_x"
  y <- data_sample$"sample_y"
  true_y <- data_sample$"true_y"
  true_CP <- data_sample$"true CPs"
  h <- 0.01
  range <- 1

  results_BIC <- main.SCHACE(
    y,
    tuning = "IC",
    IC = "BIC.Chen"
  )

  df[i] <- results_BIC$df
  lambda[i] <- results_BIC$lambda

  CP[i] <- results_BIC$`number of CPs`

  infor.result <- Informative(sample_x[results_BIC$`locations of detected CPs`], true_CP, range)
  Informative_Sets <- infor.result$"infor"
  Uninformative_Sets <- infor.result$"uninfor"

  catch_results <- catch(sample_x, Informative_Sets, true_CP, h)
  TP[i] <- catch_results$"Number_of_Catched_CPs"

  FDR[i] <- (CP[i]-TP[i])/max(CP[i], 1)
  TDR <- TP[i] / length(true_CP)

  MSE[i] <-  mean((results_BIC$`predicted y` - true_y)^2)
  
  print(i)

}

