## Install the following packages by uncmmenting and running the code
# install.packages(c("glmnet", "monomvn", "grpreg","ordinalNet","Matrix"))
# if (!require("remotes")) install.packages("remotes")
# remotes::install_github("cran/EBglmnet")
library(EBglmnet)

# Function: Penalty Regression MSE Calculator (Ridge-type)
# Implements the approach from Gertheiss & Tutz (2016) using modern packages
compute_penalty_mse <- function(X_train, y_train, true_coef, coef_list, X_test, y_test) {
  require(glmnet)
  
  # 1. Prepare group penalty weights
  group_weights <- sapply(coef_list, function(vars) {
    sqrt(length(vars))  # Default group weighting scheme
  })
  
  # 2. Create penalty factor vector aligned with X columns
  penalty_factor <- rep(1, ncol(X_train))
  names(penalty_factor) <- colnames(X_train)
  
  # Assign group weights to individual coefficients
  for(i in seq_along(coef_list)) {
    penalty_factor[coef_list[[i]]] <- group_weights[i]
  }
  
  # 3. Fit penalized model with 10-fold CV
  cv_fit <- cv.glmnet(
    X_train, y_train,
    alpha = 0,          # Ridge penalty
    penalty.factor = penalty_factor,
    standardize = FALSE # Maintain original scaling for categoricals
  )
  
  # 4. Extract coefficients at optimal lambda
  fit_coef <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
  names(fit_coef) <- colnames(X_train)
  mu <- coef(cv_fit, s = "lambda.min")[1]
  y_pred_test <- X_test %*% fit_coef + mu
  # 5. Calculate MSE metrics
  results <- list(
    # A. Prediction MSE
    prediction_mse = min(cv_fit$cvm),
    
    # B. Coefficient MSEs per categorical variable
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - fit_coef[vars])^2)
    }),
    
    # C. Full coefficient estimates
    estimates = fit_coef,
    test_mse = mean((y_test - y_pred_test) ^ 2)
  )
  
  return(results)
}

# Function: Bayesian Lasso MSE Calculator
compute_blasso_mse <- function(X_train, y_train, true_coef, coef_list, X_test, y_test) {
  require(monomvn)
  # 1. Fit Bayesian Lasso with MCMC
  fit <- blasso(X_train, y_train, 
                T = 3200,       # Total iterations (1000 post-burnin)
                thin = 1,       # No thinning
                normalize = FALSE)
  
  # 2. Manual burnin: Discard first 300 samples
  keep_samples <- 301:3200
  post_beta <- fit$beta[keep_samples, ]
  mu <- mean(fit$mu[keep_samples])
  
  # 3. Extract posterior mean coefficients
  post_mean <- colMeans(post_beta)
  names(post_mean) <- colnames(X_train)
  
  # 3.5 Checking for Convergence of MCMC
  plot(fit$beta[, 1], type = 'l', main = "MCMC Trace Plot")
  
  # 4. Calculate MSE metrics
  results <- list(
    prediction_mse = mean((y_train - (X_train %*% post_mean) - mu)^2),
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - post_mean[vars])^2)
    }),
    estimates = post_mean,
    test_mse = mean((y_test - (X_test %*% post_mean) - mu) ^ 2)
  )
  
  return(results)
}

# Function: Bayesian Elastic Net MSE Calculator
compute_ben_mse <- function(X_train, y_train, true_coef, coef_list, 
                            alpha = 0.2, lambda = 0.01) {
  require(EBglmnet)
  
  # 1. Get all variable names from design matrix
  all_vars <- colnames(X_train)
  n_vars <- length(all_vars)
  
  # 2. Fit EBglmnet model
  fit <- EBglmnet::EBglmnet(x = X_train, 
                            y = y_train,
                            family = "gaussian",
                            prior = "elastic net",
                            hyperparameters = c(alpha, lambda))
  
  # 3. Extract coefficients using predictor indices
  if(nrow(fit$fit) > 0) {
    # Get column indices from 'predictor' column (1-based)
    var_indices <- as.numeric(fit$fit[, "predictor"])
    
    # Convert indices to variable names
    nonzero_vars <- all_vars[var_indices]
    nonzero_beta <- fit$fit[, "beta"]
  } else {
    nonzero_vars <- character(0)
    nonzero_beta <- numeric(0)
  }
  
  # 4. Create full coefficient vector with zeros
  ben_coef <- setNames(rep(0, n_vars), all_vars)
  ben_coef[nonzero_vars] <- nonzero_beta
  
  # 5. Add intercept
  ben_coef <- c("(Intercept)" = fit$Intercept, ben_coef)
  
  # 6. Calculate MSE metrics
  results <- list(
    prediction_mse = mean((y_train - (X_train %*% ben_coef[-1] + ben_coef[1]))^2),
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - ben_coef[vars])^2)
    }),
    estimates = ben_coef[-1],
    test_mse = mean((y_test - (X_test %*% ben_coef[-1] + ben_coef[1])) ^ 2)
  )
  
  return(results)
}

compute_glasso_mse <- function(df_train,
                               X_train,     # your dummy‐coded training matrix
                               y_train,
                               true_coef,   # named vector of “true” coefficients, names matching colnames(X_train)
                               coef_list,   # list: each element is a character vector of column‐names belonging to one covariate
                               df_test,
                               X_test,      # your dummy‐coded test matrix (same columns / order as X_train)
                               y_test) {
  require(grpreg)
  
  # 1. Build group‐index from coef_list
  #    e.g. if coef_list = list(C1 = c("C12","C13",…),
  #                              C2 = c("C22","C23",…), …)
  #    then groups = c(1,1,…, 2,2,…, 3,3,…)
  group_ids <- rep(seq_along(coef_list), times = sapply(coef_list, length))
  
  # 2. Fit grpreg with CV
  cvfit <- cv.grpreg(X_train, y_train,
                     group   = group_ids,
                     penalty = "grLasso",
                     returnX = FALSE)
  
  # 3. Extract the coefficients at λ_min
  #    `coef(cvfit)` returns a length‐(p+1) vector: intercept + p betas
  all_coef <- coef(cvfit, lambda = cvfit$lambda.min)
  intercept <- as.numeric(all_coef[1])
  betas     <- as.numeric(all_coef[-1])
  names(betas) <- colnames(X_train)
  
  # 4. Compute per‐covariate coefficient MSE
  coef_mse <- sapply(coef_list, function(cols) {
    mean( ( true_coef[ cols ] - betas[ cols ] )^2 )
  })
  
  # 5. Cross‐validated training‐MSPE
  #    cvfit$cve is the vector of CV errors; cvfit$min is the index of λ_min
  prediction_mse <- cvfit$cve[ cvfit$min ]
  
  # 6. Compute test‐MSE
  #    ŷ_test = intercept + X_test %*% betas
  y_hat_test <- as.vector( X_test %*% betas + intercept )
  test_mse   <- mean( (y_test - y_hat_test)^2 )
  
  # 7. Return all four pieces
  return(list(
    prediction_mse  = prediction_mse,
    coefficient_mse = coef_mse,
    estimates       = betas,
    test_mse        = test_mse
  ))
}






