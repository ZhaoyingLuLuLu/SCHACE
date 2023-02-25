
CV.SCHACE <- function(A, N, sample_x, sample_y, folds, Bdf_set, methods, clambda, percent){

  folds.index = foldK(N, folds)
  size = rep(0, folds)

  for(i in 1:folds){

    size[i] = length(which(folds.index==i,arr.ind=TRUE))

  }

  cvm_1se <- rep(0, length(Bdf_set))
  cvm_min <- rep(0, length(Bdf_set))
  lambda.1se <- rep(0, length(Bdf_set))
  lambda.min <- rep(0, length(Bdf_set))

  for(k in 1:length(Bdf_set)){   #------------------ tune Bdf -------------------------

    #k = 3

    bsp_mat = splines::bs(sample_x, df = Bdf_set[k])
    X = cbind(A, bsp_mat[, 1:Bdf_set[k]])
    data = cbind(sample_y, X)

    object = glmnet::glmnet(X, sample_y, family = "gaussian", penalty.factor = c(rep(1, N-1), rep(0, Bdf_set[k])),
                    standardize = TRUE, gamma = 0, relax = TRUE)
    lambda_set = object$lambda
    nla = length(lambda_set)

    pred_y_mat_relax = matrix(, ncol = nla, nrow = length(sample_y))
    cvraw = matrix(, ncol = nla, nrow = length(sample_y))
    cvraw_foldmean = matrix(, ncol = nla, nrow = folds)

    for (i in 1:folds){  #each fold
      #i=1
      testIndexes = which(folds.index==i,arr.ind=TRUE)

      testData = data[testIndexes, ]
      testData.X = testData[, -1]
      testData.y = testData[, 1]

      trainData = data[-testIndexes, ]
      trainData.X = trainData[, -1]
      trainData.y = trainData[, 1]

      object.relaxed = glmnet::glmnet(trainData.X, trainData.y, family = "gaussian",
                              penalty.factor = c(rep(1, N-1), rep(0, Bdf_set[k])),
                              standardize = TRUE, lambda = lambda_set, gamma = 0, relax = TRUE)


      preds = predict(object.relaxed$relaxed, testData.X, s=lambda_set)
      pred_y_mat_relax[testIndexes, ] = preds

      if(methods == "MSE"){

        for (m in 1:length(lambda_set)){

          cvraw[testIndexes, m] = meansquared(testData.y, pred_y_mat_relax[testIndexes, m])
          cvraw_foldmean[i, m] = mean(cvraw[testIndexes, m])

        }

      }else if(methods == "Abs"){

        for (m in 1:length(lambda_set)){

          cvraw[testIndexes, m] = absolute(testData.y, pred_y_mat_relax[testIndexes, m])
          cvraw_foldmean[i, m] = mean(cvraw[testIndexes, m])

        }

      }else if(methods == "Trim"){

        for (m in 1:length(lambda_set)){

          cvraw_foldmean[i, m] = trimmedmse(testData.y, pred_y_mat_relax[testIndexes, m], percent)

        }

      }

    }

    cvm = apply(cvraw_foldmean, 2, weighted.mean, w = size, na.rm = TRUE)  #weighted mean

    sd_sum = apply(scale(cvraw_foldmean, cvm, FALSE)^2, 2, weighted.mean, w = size, na.rm = TRUE)

    cvsd = sqrt(sd_sum/(folds-1))

    #---- from function getOptcv.glmnet() -----
    cvmin = min(cvm, na.rm = TRUE)
    if(clambda == "lambda1se"){
      idmin = cvm <= cvmin
    }else if(clambda == "lambdamin"){
      idmin = which(abs(cvm-cvmin) <= 10^(-5))
    }

    lambda_min = max(lambda_set[idmin], na.rm = TRUE)
    idmin = match(lambda_min, lambda_set)   #index
    semin = (cvm + cvsd)[idmin]
    id1se = cvm <= semin
    #------------------------------------------

    lambda.1se[k] = max(lambda_set[id1se], na.rm = TRUE)   #lambda.1se
    lambda.min[k] = lambda_min                             #lambda.min

    #cvm_all[k] = min(cvm)
    index_lambda1se <- which(lambda_set == lambda.1se[k])
    cvm_1se[k] = cvm[index_lambda1se]
    cvm_min[k] = cvmin

  }

  if(clambda == "lambda1se"){

    optimal_lambda <- lambda.1se[which.min(cvm_1se)]
    Bdf <- Bdf_set[which.min(cvm_1se)]

  }else if(clambda == "lambdamin"){

    optimal_lambda <- lambda.min[which.min(cvm_min)]
    Bdf <- Bdf_set[which.min(cvm_min)]

  }

  results <- list("df" = Bdf, "lambda" = optimal_lambda)

  return(results)


}
