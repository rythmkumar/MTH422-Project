# Function: Penalty Regression MSE Calculator (Ridge-type)
# Implements the approach from Gertheiss & Tutz (2016) using modern packages
compute_penalty_mse <- function(X_train, y_train, true_coef, coef_list) {
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
  
  # 5. Calculate MSE metrics
  results <- list(
    # A. Prediction MSE
    prediction_mse = min(cv_fit$cvm),
    
    # B. Coefficient MSEs per categorical variable
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - fit_coef[vars])^2)
    }),
    
    # C. Full coefficient estimates
    estimates = fit_coef
  )
  
  return(results)
}

# Function: Bayesian Lasso MSE Calculator
compute_blasso_mse <- function(X_train, y_train, true_coef, coef_list) {
  require(monomvn)
  
  # 1. Fit Bayesian Lasso with MCMC
  fit <- blasso(X_train, y_train, 
                T = 1000,       # Total MCMC iterations
                burnin = 200,   # Burn-in period
                normalize = FALSE)  # Keep dummy variables unscaled
  
  # 2. Extract posterior mean coefficients (excluding intercept)
  post_mean <- colMeans(fit$beta)
  names(post_mean) <- colnames(X_train)
  
  # 2.5 Checking for Convergence of MCMC
  plot(fit$beta[, 1], type = 'l', main = "MCMC Trace Plot")
  
  # 3. Calculate MSE metrics
  results <- list(
    prediction_mse = mean((y_train - (X_train %*% post_mean))^2),
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - post_mean[vars])^2)
    }),
    estimates = post_mean
  )
  
  return(results)
}

# Function: Bayesian Elastic Net MSE Calculator
compute_ben_mse <- function(X_train, y_train, true_coef, coef_list, 
                            alpha = 0.5, # Elastic Net mixing parameter
                            n_iter = 1000) {
  require(BGLR)
  
  # 1. Prepare data for BGLR (requires centered response)
  y_centered <- y_train - mean(y_train)
  
  # 2. Set up Elastic Net prior
  ETA <- list(EN = list(X = X_train, 
                        model = "EN",
                        lambda = 0.1,   # Initial regularization
                        alpha = alpha)) # α = 0 (ridge), 1 (lasso)
  
  # 3. Fit model
  fit <- BGLR(y = y_centered,
              ETA = ETA,
              nIter = n_iter,
              burnIn = 200,
              verbose = FALSE)
  
  # 4. Extract coefficients (re-center intercept)
  ben_coef <- c(fit$mu + mean(y_train), fit$ETA$EN$b)
  names(ben_coef) <- c("(Intercept)", colnames(X_train))
  
  # 5. Calculate MSE metrics
  results <- list(
    prediction_mse = mean((y_train - (X_train %*% ben_coef[-1] + ben_coef[1]))^2),
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - ben_coef[vars])^2)
    }),
    estimates = ben_coef[-1]
  )
  
  return(results)
}

compute_glasso_mse <- function(df_train, y_train, true_coef, coef_list) {
  require(grpreg)
  
  # 1. Preprocessing: Apply polynomial contrasts to ordinal variables (C1-C4)
  df_processed <- df_train
  ordinal_predictors <- c("C1", "C2", "C3", "C4")
  
  # Apply polynomial contrasts to ordinal variables
  for(pred in ordinal_predictors) {
    n_levels <- length(levels(df_processed[[pred]]))
    contrasts(df_processed[[pred]]) <- contr.poly(n_levels)
  }
  
  # 2. Create design matrix with updated contrasts
  X_glasso <- model.matrix(~C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8, 
                           data = df_processed)[,-1]
  
  # 3. Define group structure (1 group per predictor)
  groups <- rep(1:8, times = sapply(coef_list, length))
  
  # 4. Fit Group Lasso with cross-validation
  cv_fit <- cv.grpreg(X_glasso, y_train, 
                      group = groups,
                      penalty = "grLasso",
                      returnX = FALSE)
  
  # 5. Extract coefficients at optimal lambda
  fit_coef <- coef(cv_fit, lambda = cv_fit$lambda.min)
  fit_coef <- fit_coef[-1]  # Remove intercept
  
  # 6. Map coefficients back to original dummy names
  orig_names <- colnames(model.matrix(~., df_train))[-1]
  aligned_coef <- setNames(numeric(length(orig_names)), orig_names)
  aligned_coef[names(fit_coef)] <- fit_coef
  
  # 7. Calculate MSE metrics
  results <- list(
    prediction_mse = cv_fit$cve[cv_fit$min],
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - aligned_coef[vars])^2)
    }),
    estimates = aligned_coef
  )
  
  return(results)
}

compute_glap_mse <- function(df_train, y_train, true_coef, coef_list) {
  require(ordinalNet)
  require(Matrix)
  
  # 1. Preprocessing: Convert ordinal predictors to ordered factors
  df_processed <- df_train
  ordinal_preds <- c("C1", "C2", "C3", "C4")
  df_processed[ordinal_preds] <- lapply(df_processed[ordinal_preds], ordered)
  
  # 2. Create design matrix with ordinal-aware encoding
  X_glap <- model.matrix(~C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8, 
                         data = df_processed)[,-1]
  
  # 3. Create penalty matrices for ordinal variables
  ord_mats <- lapply(ordinal_preds, function(pred) {
    n_levels <- nlevels(df_processed[[pred]])
    get_laplacian(n_levels - 1)  # From ordinalNet
  })
  
  # 4. Define penalty factors for grouped Laplacian
  penalty_mats <- do.call("bdiag", ord_mats)
  groups <- c(rep(1:4, each = 1),  # Ordinal groups with Laplacian
              rep(5:8, times = lengths(coef_list[5:8])))  # Nominal groups
  
  # 5. Fit ordinal-aware elastic net with Laplacian penalty
  fit <- ordinalNet(X_glap, y_train,
                    family = "gaussian",
                    link = "identity",
                    alpha = 0.2,  # Balance L1/L2
                    lambdaValues = NULL,
                    penaltyFactors = c(rep(1, 4), rep(0, 4)),  # Penalize ordinals
                    standardize = FALSE)
  
  # 6. Extract coefficients (excluding intercept)
  fit_coef <- coef(fit)[-1]
  names(fit_coef) <- colnames(X_glap)
  
  # 7. Calculate MSE metrics
  results <- list(
    prediction_mse = mean((y_train - (X_glap %*% fit_coef))^2),
    coefficient_mse = sapply(coef_list, function(vars) {
      mean((true_coef[vars] - fit_coef[vars])^2)
    }),
    estimates = fit_coef
  )
  
  return(results)
}

# Helper function to create Laplacian matrix
get_laplacian <- function(n) {
  D <- diff(diag(n), differences = 1)
  t(D) %*% D
}


