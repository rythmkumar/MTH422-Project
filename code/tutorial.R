# Install and load necessary packages
# install.packages(c("bayesm","mcclust"))
# install.packages("remotes")
# remotes::install_github("cran/GreedyEPL")
# install.packages("https://cran.r-project.org/src/contrib/Archive/effectFusion/effectFusion_1.1.3.tar.gz",
#                  repos = NULL, type = "source")

# Set seed for reproducibility
set.seed(10) ## please give us A*

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


# n_datasets <- 1 ## change accordingly 
# for (i in 1:n_datasets) {
#   data_sim <- generate_data()
#   y_train <- data_sim$y
#   df_train <- data_sim$df
#   X_train <- data_sim$X
#   ## PLease put the models code here for getting the MSE's values
#   ben_mse <- compute_ben_mse(X_train, y_train, true_coef, coef_list)
#   blasso_mse <- compute_blasso_mse(X_train, y_train, true_coef, coef_list)
#   full_mse <- compute_full_mse(df_train, X_train, y_train, true_coef, coef_list)
#   fusion_mse <- compute_fusion_mse(df_train, X_train, y_train, true_coef, coef_list)
#   # glap_mse <- compute_glap_mse(df_train, y_train, true_coef, coef_list)
#   glasso_mse <- compute_glasso_mse(df_train, y_train, true_coef, coef_list)
#   penalty_mse <- compute_penalty_mse(X_train, y_train, true_coef, coef_list)
#   true_mse <- compute_true_mse(df_train, X_train, y_train, true_coef, coef_list)
# 
# }

# priorr : mu * sigma.2 * beta (with hyperparameters: gamma^2 * Q^-1(delta) * (c / 2))
# Q^-1 is the prior precision matrix
# mu follows flat normal distribution N (0, M0)
# sigma.2 follows inverse gamma (gh0, Gh0)
# delta follows det(q)^-1/2 *root  r ^{sum (1 - delta_{kj})}
# 
# mu follows flat normal distribution N (0, M0)
# mu follows flat normal distribution N (0, M0)



### Generating plots in figure 4 and 5

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


# For coefficient MSE
coef_data <- data.frame(Method = character(), Covariate = integer(), MSE = numeric())
for (m in methods) {
  # print(m)
  for (c in 1:n_covariates) {
    # print(c)
    temp <- data.frame(Method = m,
                       Covariate = c,
                       MSE = coef_mse_list[[m]][[c]])
    rownames(temp) <- NULL
    coef_data <- rbind(coef_data, temp)
  }
}

# coef_data

# For prediction MSE
prediction_data <- data.frame(Method = character(), MSE = numeric())
for (m in methods) {
  temp <- data.frame(Method = m, MSE = train_mse_list[[m]])
  rownames(temp) <- NULL
  prediction_data <- rbind(prediction_data, temp)
}

# pred_data

# Close all open graphics devices
while (!is.null(dev.list())) dev.off()

# Load required libraries
library(ggplot2)
library(dplyr)

coef_data$Covariate <- as.factor(coef_data$Covariate)

# Loop through each covariate
for (i in 1:8) {
  p <- ggplot(coef_data %>% filter(Covariate == i), aes(x = Method, y = MSE)) +
    geom_boxplot(fill = "#56B4E9", color = "black", outlier.shape = 1, outlier.alpha = 0.5) +
    labs(title = paste("Covariate", i),
         x = "Method",
         y = "Coefficient MSE") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.y = element_line(color = "gray80")
    )
  
  print(p)
  # Optional: Pause between plots if needed
  readline(prompt = "Press [Enter] to continue to next plot...")
}


prediction_data$Method <- factor(prediction_data$Method,
                                 levels = methods)

# Boxplot for prediction MSE
ggplot(prediction_data, aes(x = Method, y = MSE)) +
  geom_boxplot(outlier.shape = 1,  # hollow circles
               outlier.size = 2,
               fill = "white",
               color = "black",
               width = 0.6,
               fatten = 1) +
  labs(title = "Simulation Study: MSPE of Training Data",
       x = NULL,
       y = "Mean Squared Prediction Error (MSPE)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray85")
  )

# save(train_beta_list, file = "training_betas.RData")
# save(train_mse_list, file = "training_mse.RData")
# save(data_sim_test, file = "test_data.RData")
# save(coef_mse_list, file = "coef_mse_list.RData")

## Boxplot for test data
# For prediction MSE

test_data <- data.frame(Method = character(), MSE = numeric())
for (m in methods) {
  temp <- data.frame(Method = m, MSE = test_mse_list[[m]])
  rownames(temp) <- NULL
  test_data <- rbind(test_data, temp)
}
# Boxplot for prediction MSE
ggplot(test_data, aes(x = Method, y = MSE)) +
  geom_boxplot(outlier.shape = 1,  # hollow circles
               outlier.size = 2,
               fill = "white",
               color = "black",
               width = 0.6,
               fatten = 1) +
  labs(title = "Simulation Study: MSPE of Testing Data",
       x = NULL,
       y = "Mean Squared Prediction Error (MSPE)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "gray85")
  )
