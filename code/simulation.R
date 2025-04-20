library(EBglmnet)

# Set seed for reproducibility
set.seed(10) ## please give us A*

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

# Function for Fusion Model evaluation
compute_fusion_mse <- function(df_train, X_train, y_train, true_coef, coef_list, df_test, X_test, y_test) {
  require(effectFusion)
  
  # 1. Fit fusion model with spike and slab prior
  fusion_fit <- effectFusion(
    y = y_train,
    X = df_train, 
    types = sim1$types,
    method = "SpikeSlab",
    prior = list(r = 20000, G0 = 20),
    mcmc = list(M = 10000, burnin = 5000),
    modelSelection = "binder"
  )
  
  # 2. Extract posterior means of coefficients
  invisible_output <- capture.output(coef_det <- summary(fusion_fit)[2:nrow(summary(fusion_fit)),])
  fit_coef <- coef_det[,1]
  invisible_output <- capture.output(mu_coef <- summary(fusion_fit)[1,1])
  names(fit_coef)<- names(true_coef)
  
  # 3. Calculate prediction MSE
  y_pred <- mu_coef + as.numeric(X_train %*% fit_coef)
  prediction_mse <- mean((y_train - y_pred)^2)
  
  y_test_pred <- mu_coef + as.numeric(X_test %*% fit_coef)
  test_mse <- mean((y_test - y_test_pred) ^ 2)
  
  # 4. Calculate coefficient MSEs per categorical variable
  coefficient_mse <- sapply(coef_list, function(vars) {
    mean((true_coef[vars] - fit_coef[vars])^2)
  })
  
  # 5. Return results
  results <- list(
    prediction_mse = prediction_mse,
    coefficient_mse = coefficient_mse,
    estimates = fit_coef,
    mu = mu_coef,
    test_mse = test_mse
  )
  
  return(results)
}

compute_full_mse <- function(df_train, X_train, y_train, true_coef, coef_list, df_test, X_test, y_test) {
  require(effectFusion)
  
  # 1. Fit full model (no fusion)
  full_fit <- effectFusion(
    y = y_train,
    X = df_train, 
    types = sim1$types,
    method = NULL,  # No fusion
    mcmc = list(M = 3000, burnin = 1000)  # Using shorter chain for full model
  )
  
  # 2. Extract posterior means of coefficients
  invisible_output <- capture.output(coef_det <- summary(full_fit)[2:nrow(summary(full_fit)),])
  fit_coef <- coef_det[,1]
  invisible_output <- capture.output(mu_coef <- summary(full_fit)[1,1])
  names(fit_coef) <- names(true_coef)
  # print("full")
  # print(fit_coef)
  # 3. Calculate prediction MSE
  y_pred <- mu_coef + as.numeric(X_train %*% fit_coef)
  prediction_mse <- mean((y_train - y_pred)^2)
  
  y_pred_test <- mu_coef + as.numeric(X_test %*% fit_coef)
  test_mse <- mean((y_test - y_pred_test) ^ 2)
  
  # 4. Calculate coefficient MSEs per categorical variable
  coefficient_mse <- sapply(coef_list, function(vars) {
    mean((true_coef[vars] - fit_coef[vars])^2)
  })
  
  # 5. Return results
  results <- list(
    prediction_mse = prediction_mse,
    coefficient_mse = coefficient_mse,
    estimates = fit_coef,
    mu = mu_coef,
    test_mse = test_mse
  )
  
  return(results)
}

generate_X_true <- function(df_train, true_coef) {
  # Create a copy of the original data frame
  X_true <- df_train
  
  # C1: Fusion pattern (0,1,1,2,2,4,4)
  # Level 1 = reference/0, Level 2 = 0 (should be fused with level 1), 
  # Levels 3-4 = 1, Levels 5-6 = 2, Levels 7-8 = 4
  X_true$C1 <- factor(ifelse(df_train$C1 %in% c(1, 2), 1,  # Fuse levels 1-2 (both have effect 0)
                             ifelse(df_train$C1 %in% c(3, 4), 3,    # Fuse levels 3-4 (effect = 1)
                                    ifelse(df_train$C1 %in% c(5, 6), 5,    # Fuse levels 5-6 (effect = 2)
                                           ifelse(df_train$C1 %in% c(7, 8), 7,    # Fuse levels 7-8 (effect = 4)
                                                  df_train$C1)))))
  
  # C2: All coefficients are 0, fuse all levels
  X_true$C2 <- factor(rep(1, nrow(df_train)))
  
  # C3: Fusion pattern (0,0,-2,-2)
  # Levels 1-2 = 0 (should be fused), Levels 3-4 = -2
  X_true$C3 <- factor(ifelse(df_train$C3 %in% c(1, 2), 1,  # Fuse levels 1-2 (both have effect 0)
                             ifelse(df_train$C3 %in% c(3, 4), 3,    # Fuse levels 3-4 (effect = -2)
                                    df_train$C3)))
  
  # C4: All coefficients are 0, fuse all levels
  X_true$C4 <- factor(rep(1, nrow(df_train)))
  
  # C5: Fusion pattern (0,0,1,1,1,1,-2,-2)
  # Levels 1-2 = 0 (should be fused), Levels 3-6 = 1, Levels 7-8 = -2
  X_true$C5 <- factor(ifelse(df_train$C5 %in% c(1, 2), 1,  # Fuse levels 1-2 (both have effect 0)
                             ifelse(df_train$C5 %in% c(3, 4, 5, 6), 3,  # Fuse levels 3-6 (effect = 1)
                                    ifelse(df_train$C5 %in% c(7, 8), 7,    # Fuse levels 7-8 (effect = -2)
                                           df_train$C5))))
  
  # C6: All coefficients are 0, fuse all levels
  X_true$C6 <- factor(rep(1, nrow(df_train)))
  
  # C7: Fusion pattern (0,0,2,2)
  # Levels 1-2 = 0 (should be fused), Levels 3-4 = 2
  X_true$C7 <- factor(ifelse(df_train$C7 %in% c(1, 2), 1,  # Fuse levels 1-2 (both have effect 0)
                             ifelse(df_train$C7 %in% c(3, 4), 3,    # Fuse levels 3-4 (effect = 2)
                                    df_train$C7)))
  
  # C8: All coefficients are 0, fuse all levels
  X_true$C8 <- factor(rep(1, nrow(df_train)))
  
  # Ensure all columns are factors
  X_true <- as.data.frame(lapply(X_true, factor))
  
  return(X_true)
}

compute_true_mse <- function(df_train, X_train, y_train, true_coef, coef_list, X_test, y_test, df_test) {
  require(effectFusion)
  
  # 1. Generate X_true using fusion patterns
  X_true <- generate_X_true(df_train, true_coef)
  X_test_true <- generate_X_true(df_test, true_coef)
  
  # 2. Fit model with the true fusion structure
  true_fit <- effectFusion(
    y = y_train, 
    X = X_true,
    types = sim1$types,
    method = NULL,
    mcmc = list(M = 3000, burnin = 1000)
  )
  
  # 3. Extract posterior means of coefficients
  invisible_output <- capture.output(coef_det <- summary(true_fit)[2:nrow(summary(true_fit)),])
  fit_coef_true <- coef_det[,1]
  invisible_output <- capture.output(mu_coef <- summary(true_fit)[1,1])
  
  # 4. Create design matrix from X_true for prediction (including only variables with multiple levels)
  X_true_matrix <- model.matrix(~ C1 + C3 + C5 + C7, data = X_true)[, -1]
  X_test_true_matrix <- model.matrix(~ C1 + C3 + C5 + C7, data = X_test_true)[, -1]
  
  # 5. Calculate prediction MSE using refitted coefficients
  names(fit_coef_true) <- colnames(X_true_matrix)
  y_pred_true <- mu_coef + as.numeric(X_true_matrix %*% fit_coef_true)
  prediction_mse <- mean((y_train - y_pred_true)^2)
  
  y_pred_true_test <- mu_coef + as.numeric(X_test_true_matrix %*% fit_coef_true)
  test_mse <- mean((y_test - y_pred_true_test)^2)
  
  # print("true")
  # print("")
  # 6. Map estimated coefficients back to original coefficient space
  # Initialize with zeros
  fit_coef <- rep(0, length(true_coef))
  names(fit_coef) <- names(true_coef)
  
  # Map C1 coefficients (based on fusion pattern (0,1,1,2,2,4,4))
  # Reference level (C12) stays at 0
  if("C13" %in% names(fit_coef)) fit_coef["C13"] <- fit_coef_true["C12"] # Level 2
  if("C14" %in% names(fit_coef)) fit_coef["C14"] <- fit_coef_true["C13"] # Level 3-4 fused
  if("C13" %in% names(fit_coef)) fit_coef["C13"] <- fit_coef_true["C13"] # Level 3-4 fused
  if("C15" %in% names(fit_coef)) fit_coef["C15"] <- fit_coef_true["C15"] # Level 5-6 fused
  if("C16" %in% names(fit_coef)) fit_coef["C16"] <- fit_coef_true["C15"] # Level 5-6 fused
  if("C17" %in% names(fit_coef)) fit_coef["C17"] <- fit_coef_true["C17"] # Level 7-8 fused
  if("C18" %in% names(fit_coef)) fit_coef["C18"] <- fit_coef_true["C17"] # Level 7-8 fused
  
  # C2 coefficients all remain 0 (all fused to reference)
  
  # Map C3 coefficients (based on fusion pattern (0,0,-2,-2))
  # Reference level (C32) stays at 0
  # C33 and C34 are fused and share coefficient
  if("C33" %in% names(fit_coef)) fit_coef["C33"] <- fit_coef_true["C33"] 
  if("C34" %in% names(fit_coef)) fit_coef["C34"] <- fit_coef_true["C33"]
  
  # C4 coefficients all remain 0 (all fused to reference)
  
  # Map C5 coefficients (based on fusion pattern (0,0,1,1,1,1,-2,-2))
  # Reference level (C52) stays at 0
  # C53-C56 are fused and share coefficient
  if("C53" %in% names(fit_coef)) fit_coef["C53"] <- fit_coef_true["C53"]
  if("C54" %in% names(fit_coef)) fit_coef["C54"] <- fit_coef_true["C53"]
  if("C55" %in% names(fit_coef)) fit_coef["C55"] <- fit_coef_true["C53"]
  if("C56" %in% names(fit_coef)) fit_coef["C56"] <- fit_coef_true["C53"]
  # C57-C58 are fused and share coefficient
  if("C57" %in% names(fit_coef)) fit_coef["C57"] <- fit_coef_true["C57"]
  if("C58" %in% names(fit_coef)) fit_coef["C58"] <- fit_coef_true["C57"]
  
  # C6 coefficients all remain 0 (all fused to reference)
  
  # Map C7 coefficients (based on fusion pattern (0,0,2,2))
  # Reference level (C72) stays at 0
  # C73 and C74 are fused and share coefficient
  if("C73" %in% names(fit_coef)) fit_coef["C73"] <- fit_coef_true["C73"]
  if("C74" %in% names(fit_coef)) fit_coef["C74"] <- fit_coef_true["C73"]
  
  # C8 coefficients all remain 0 (all fused to reference)
  
  # 7. Calculate coefficient MSEs per categorical variable
  coefficient_mse <- sapply(coef_list, function(vars) {
    mean((true_coef[vars] - fit_coef[vars])^2)
  })
  
  # 8. Return results
  results <- list(
    prediction_mse = prediction_mse,
    coefficient_mse = coefficient_mse,
    estimates = fit_coef,
    X_true = X_true,
    test_mse = test_mse
  )
  # print(fit_coef)
  return(results)
}

# ------------------- Simulation Study -----------------------

# Define parameters
n_obs <- 500
mu <- 1
sigma <- 1

# Define probabilities for levels
pi_8 <- c(0.1, 0.1, 0.2, 0.05, 0.2, 0.1, 0.2, 0.05)
pi_4 <- c(0.1, 0.4, 0.2, 0.3)

# Define true coefficients for each predictor's dummy variables (excluding intercept)
# first 4 (1, 2, 3, 4) are ordinal while the remaining 4 are nominal
true_coef <- c(
  "C12" = 0, "C13" = 1, "C14" = 1, "C15" = 2, "C16" = 2, "C17" = 4, "C18" = 4,
  "C22" = 0, "C23" = 0, "C24" = 0, "C25" = 0, "C26" = 0, "C27" = 0, "C28" = 0,
  "C32" = 0, "C33" = -2, "C34" = -2,
  "C42" = 0, "C43" = 0, "C44" = 0,
  "C52" = 0, "C53" = 1, "C54" = 1, "C55" = 1, "C56" = 1, "C57" = -2, "C58" = -2,
  "C62" = 0, "C63" = 0, "C64" = 0, "C65" = 0, "C66" = 0, "C67" = 0, "C68" = 0,
  "C72" = 0, "C73" = 2, "C74" = 2,
  "C82" = 0, "C83" = 0, "C84" = 0
)

# Define list of coefficients for each predictor
coef_list <- list(
  "C1" = c("C12", "C13", "C14", "C15", "C16", "C17", "C18"),
  "C2" = c("C22", "C23", "C24", "C25", "C26", "C27", "C28"),
  "C3" = c("C32", "C33", "C34"),
  "C4" = c("C42", "C43", "C44"),
  "C5" = c("C52", "C53", "C54", "C55", "C56", "C57", "C58"),
  "C6" = c("C62", "C63", "C64", "C65", "C66", "C67", "C68"),
  "C7" = c("C72", "C73", "C74"),
  "C8" = c("C82", "C83", "C84")
)

# Function to generate a single dataset
generate_data <- function() {
  C1 <- factor(sample(1:8, n_obs, TRUE, pi_8), levels = 1:8)
  C2 <- factor(sample(1:8, n_obs, TRUE, pi_8), levels = 1:8)
  C3 <- factor(sample(1:4, n_obs, TRUE, pi_4), levels = 1:4)
  C4 <- factor(sample(1:4, n_obs, TRUE, pi_4), levels = 1:4)
  C5 <- factor(sample(1:8, n_obs, TRUE, pi_8), levels = 1:8)
  C6 <- factor(sample(1:8, n_obs, TRUE, pi_8), levels = 1:8)
  C7 <- factor(sample(1:4, n_obs, TRUE, pi_4), levels = 1:4)
  C8 <- factor(sample(1:4, n_obs, TRUE, pi_4), levels = 1:4)
  
  df <- data.frame(C1, C2, C3, C4, C5, C6, C7, C8)
  X <- model.matrix(~C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8, df)  # No intercept in X
  X <- X[, -1]
  epsilon <- rnorm(n_obs, 0, sigma)
  y <- mu + as.vector(X %*% true_coef) + epsilon  # Intercept added separately
  return(list(y = y, df = df,X= X))
}


# prior : mu * sigma.2 * beta (with hyperparameters: gamma^2 * Q^-1(delta) * (c / 2))
# Q^-1 is the prior precision matrix
# mu follows flat normal distribution N (0, M0)
# sigma.2 follows inverse gamma (gh0, Gh0)
# delta follows det(q)^-1/2 *root  r ^{sum (1 - delta_{kj})}
# 
# mu follows flat normal distribution N (0, M0)
# mu follows flat normal distribution N (0, M0)


# Assuming methods list
methods <- c("Full", "Fusion", "Penalty", "BLasso", "BEN", "GLasso", "True") # ADD SGL, BSGS
n_methods <- length(methods)
n_covariates <- 8
n_betas <- 40

# Create a list to store coefficient MSEs per method and per covariate
coef_mse_list <- vector("list", length = n_methods)
train_beta_list <- vector("list", length = n_methods)

names(coef_mse_list) <- methods
for (m in methods) {
  coef_mse_list[[m]] <- vector("list", length = n_covariates)
  train_beta_list[[m]] <- vector("list", length = n_betas)
  for (c in 1:n_covariates) {
    # print(c)
    coef_mse_list[[m]][[c]] <- numeric(length = 0)
    # print(coef_mse_list[[m]][[c]])
  }
  for (cii in 1:n_betas) {
    train_beta_list[[m]][[cii]] <- numeric(length = 0)
  }
}

# Create a list to store prediction MSEs
train_mse_list <- list()
for (m in methods) {
  train_mse_list[[m]] <- numeric()
}

n_datasets <- 100

## testing for new dataset
data_sim_test <- generate_data()
y_test <- data_sim_test$y
df_test <- data_sim_test$df
X_test <- data_sim_test$X

test_mse_list <- list()
for (m in methods) {
  test_mse_list[[m]] <- numeric()
}

# train_data <- vector("list", length = n_datasets)

for (i in 1:n_datasets) {
  print(paste("dataset", i))
  data_sim <- generate_data()
  y_train <- data_sim$y
  df_train <- data_sim$df
  X_train <- data_sim$X
  
  # Assuming true_coef and coef_list are defined
  mse_results <- list(
    Full = compute_full_mse(df_train, X_train, y_train, true_coef, coef_list, df_test, X_test, y_test),
    Fusion = compute_fusion_mse(df_train, X_train, y_train, true_coef, coef_list, df_test, X_test, y_test),
    Penalty = compute_penalty_mse(X_train, y_train, true_coef, coef_list, X_test, y_test),
    BLasso = compute_blasso_mse(X_train, y_train, true_coef, coef_list, X_test, y_test),
    BEN = compute_ben_mse(X_train, y_train, true_coef, coef_list, X_test, y_test),
    # GLasso = compute_glasso_mse(df_train, y_train, true_coef, coef_list, df_test, X_test, y_test),
    GLasso = compute_glasso_mse(df_train, X_train, y_train, true_coef, coef_list, df_test, X_test, y_test),
    # GLap = compute_glap_mse(df_train, y_train, true_coef, coef_list), # Uncomment if needed
    # SGL = compute_sgl_mse(df_train, y_train, true_coef, coef_list),
    # BSGS = compute_bsgs_mse(df_train, y_train, true_coef, coef_list),
    True = compute_true_mse(df_train, X_train, y_train, true_coef, coef_list, X_test, y_test, df_test)
    
  )
  # print(mse_results$Full)
  # Store values
  # train_data[i] <- data_sim
  for (m in names(mse_results)) {
    train_mse_list[[m]] <- c(train_mse_list[[m]], mse_results[[m]]$prediction_mse)
    
    for (c in 1:n_covariates) {
      coef_mse_list[[m]][[c]] <- c(coef_mse_list[[m]][[c]], mse_results[[m]]$coefficient_mse[c])
    }
    for (cii in 1:n_betas) {
      train_beta_list[[m]][[cii]] <- c(train_beta_list[[m]][[cii]], mse_results[[m]]$estimates[cii])
    }
    test_mse_list[[m]] <- c(test_mse_list[[m]], mse_results[[m]]$test_mse)
  }
  
}

save(train_beta_list, file = "training_betas.RData")
save(train_mse_list, file = "training_mse.RData")
save(data_sim_test, file = "test_data.RData")
save(coef_mse_list, file = "coef_mse_list.RData")
save(test_mse_list, file =  "test_mse.RData")

