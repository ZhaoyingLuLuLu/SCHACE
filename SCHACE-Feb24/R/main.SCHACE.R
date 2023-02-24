
main.SCHACE <- function(y, folds = 3, Bdf_set=seq(3, 15, 1),
                     lambda = NULL, Bdf = NULL,
                     tuning = c("IC", "crossvalidation"), IC = c("BIC", "BIC.Chen"),
                     methods = c("MSE", "Abs", "Trim"),
                     clambda = "lambdamin", percent = 0.1){

  sample_x = 1:length(y)
  N = length(y)

  A = matrix(rep(1, N*N), ncol = N)
  A_whole = A*lower.tri(A)
  A = A_whole[, 1:(N-1)]

  if (is.null(lambda) & is.null(Bdf)){

    fit = switch(tuning,
                 "IC" = IC.SCHACE(A, N, sample_x, y, Bdf_set, IC),
                 "crossvalidation" = CV.SCHACE(A, N, sample_x, y, folds, Bdf_set, methods, clambda, percent)
    )

    Bdf <- fit$df
    lambda <- fit$lambda

  }else if(is.null(lambda) & !is.null(Bdf)){

    stop("Please specify lambda as well.")

  }else if(!is.null(lambda) & is.null(Bdf)){

    stop("Please specify B-splines Degree of Freedom as well.")

  }

  bsp_mat  <-  bs(sample_x, df = Bdf)
  X <- cbind(A, bsp_mat[, 1:Bdf])

  object = glmnet(X, y, family = "gaussian", penalty.factor = c(rep(1, N-1), rep(0, Bdf)),
                  standardize = TRUE, gamma = 0, relax = TRUE)

  beta = as.vector(coef(object, s = lambda, gamma = 0, relax = TRUE))

  datapoints = beta[2:N]
  location = which(datapoints != 0)
  changepoints = sample_x[location]
  number_CPs = length(changepoints)
  jumpsize <- datapoints[location]

  pred_y <- predict(object$relaxed, newx = X, type = "response", s = lambda)

  results <- list("df" = Bdf, "lambda" = lambda,
                  "predicted y" = pred_y, "beta" = beta,
                  "number of CPs" = number_CPs, "locations of detected CPs" = changepoints, "jump size" = jumpsize)

  return(results)

}
