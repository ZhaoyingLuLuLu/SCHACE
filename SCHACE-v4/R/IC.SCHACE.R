
IC.SCHACE <- function(A, N, sample_x, sample_y, Bdf_set, IC.method){

  nla = 100
  lambda_matrix = matrix(, nrow = length(Bdf_set), ncol = nla)

  BIC = matrix(, nrow = length(Bdf_set), ncol = nla)
  trainerr = matrix(, nrow = length(Bdf_set), ncol = nla)
  M_hat = matrix(, nrow = length(Bdf_set), ncol = nla)

  for (j in 1:length(Bdf_set)){

    #j = 1

    bsp_mat = splines::bs(sample_x, df = Bdf_set[j])

    X = cbind(A, bsp_mat[, 1:Bdf_set[j]])

    object = glmnet::glmnet(X, sample_y, family = "gaussian", penalty.factor = c(rep(1, N-1), rep(0, Bdf_set[j])),
                    standardize = TRUE, gamma = 0, relax = TRUE)
    lambda_set = object$lambda
    lambda_matrix[j, ] = lambda_set

    coefficient = as.matrix(coef(object$relaxed, s = lambda_set))

    pred_y = predict(object$relaxed, newx = X, type = "response", s = lambda_set)


    for(i in 1:length(lambda_set)){

      M_hat[j, i] = length( which(abs(coefficient[, i]) > 1e-8 ) )
      trainerr[j, i] = mean((pred_y[, i] - sample_y)^2)

      BIC[j, i] = N*log(trainerr[j, i]) + log(N)*M_hat[j, i]


    }

  }

  if(IC.method == "BIC.Chen"){

    BIC = BIC + 2*M_hat*log(ncol(X))

  }

  xy = as.matrix(which(BIC == min(BIC), arr.ind=T))

  Bdf = Bdf_set[xy[1, 1]]
  optimal_lambda = lambda_matrix[xy[1, 1], xy[1, 2]]      # the largest lambda with the same BIC


  results <- list("df" = Bdf, "lambda" = optimal_lambda)


  return(results)

}
