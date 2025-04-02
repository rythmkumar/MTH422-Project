# Install and load necessary packages
# install.packages("rstanarm")
# library(rstanarm)
install.packages(c("rstanarm", "monomvn", "EBglmnet", "grpreg"))
library(rstanarm)
library(monomvn)
library(EBglmnet)
library(grpreg)
library(rjags)

# Set seed for reproducibility
set.seed(123)

# Define parameters
# n_datasets <- 100
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
  
  return(list(y = y, df = df, X = X))
}

fit_glasso <- function(y, X) {
  # data_sim <- generate_data()
  # y <- data_sim$y
  # df_train <- data_sim$df
  # X <- data_sim$X
  # 
  group <- rep(1:8, times = sapply(coef_list, length))
  fit <- grpreg(X, y, group, penalty = "grLasso", nlambda = 100)
  intercept <- fit$beta[, 100][1]
  coefs <- fit$beta[, 100][-1]  # Only dummy coefficients
  names(coefs) <- colnames(X)
  return(list(coefs = c("(Intercept)" = intercept, coefs), 
              pred = function(X_new) intercept + X_new %*% coefs))
}


mse <- 0
n_datasets <- 1
for (i in 1:n_datasets) {
  data_sim <- generate_data()
  y_train <- data_sim$y
  df_train <- data_sim$df
  X_train <- data_sim$X
  
  
  ## glasso
  # fit_res <- fit_glasso(y_train, X_train)
  # coefs <- fit_res$coefs
  
}

for (h in seq_along(coef_list)) {
  coefs_h <- coef_list[[h]]
  print(coefs_h)
  mse <- mse + mean((coefs[coefs_h] - true_coef[coefs_h]) ^ 2)
}
pr
# priorr : mu * sigma.2 * beta (with hyperparameters: gamma^2 * Q^-1(delta) * (c / 2))
# Q^-1 is the prior precision matrix
# mu follows flat normal distribution N (0, M0)
# sigma.2 follows inverse gamma (gh0, Gh0)
# delta follows det(q)^-1/2 *root  r ^{sum (1 - delta_{kj})}
# 
# mu follows flat normal distribution N (0, M0)
# mu follows flat normal distribution N (0, M0)
