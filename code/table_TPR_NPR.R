library(dplyr)

load("training_betas.RData")

original_levels <- c(7, 7, 3, 3, 7, 7, 3, 3)
n_covs <- length(original_levels)
n_datasets <- 100
threshold <- 1e-6

# Compute start and end indices for original covariate levels
start_idx <- cumsum(c(1, original_levels[-length(original_levels)]))
end_idx <- cumsum(original_levels)

# Construct new covariate index list (48 levels total)
covariate_indices <- list()
new_index <- 1
for (i in seq_along(original_levels)) {
  covariate_indices[[i]] <- new_index:(new_index + original_levels[i])
  new_index <- new_index + original_levels[i] + 1
}

# Add zero-valued level to each covariate in every model
pad_beta_list <- function(beta_list, method_name) {
  extended_list <- vector("list", sum(original_levels) + length(original_levels))  # 40 + 8 = 48
  read_index <- 1
  write_index <- 1
  for (i in seq_along(original_levels)) {
    for (j in 1:original_levels[i]) {
      extended_list[[write_index]] <- beta_list[[read_index]]
      read_index <- read_index + 1
      write_index <- write_index + 1
    }
    # Add extra (zero or NA) level
    if (method_name == "True") {
      extended_list[[write_index]] <- rep(0, n_datasets)
    } else {
      extended_list[[write_index]] <- rep(0, n_datasets)
    }
    write_index <- write_index + 1
  }
  return(extended_list)
}


# Get pairwise difference vector (0 = same, 1 = different)
get_pairwise_diff <- function(beta_vec, cov_indices) {
  diffs <- c()
  for (indices in cov_indices) {
    pairs <- combn(indices, 2)
    for (i in 1:ncol(pairs)) {
      diff_val <- abs(beta_vec[pairs[1, i]] - beta_vec[pairs[2, i]])
      diffs <- c(diffs, as.numeric(diff_val > threshold))
    }
  }
  return(diffs)
}

# Metrics calculation
get_metrics <- function(truth, pred) {
  TP <- sum(truth == 1 & pred == 1)
  TN <- sum(truth == 0 & pred == 0)
  FP <- sum(truth == 0 & pred == 1)
  FN <- sum(truth == 1 & pred == 0)
  TPR <- if ( !is.na(TP) && !is.na(FN) && ((TP + FN) > 0)) TP / (TP + FN) else NA
  TNR <- if (!is.na(TN) && !is.na(FP) && ((TN + FP) > 0)) TN / (TN + FP) else NA
  PPV <- if (!is.na(TP) && !is.na(FP)&& ((TP + FP) > 0)) TP / (TP + FP) else NA
  NPV <- if ( !is.na(TN) && !is.na(FN) && ((TN + FN) > 0)) TN / (TN + FN) else NA
  return(c(TPR = TPR, TNR = TNR, PPV = PPV, NPV = NPV))
}

# Extend all beta lists with one more level per covariate
methods <- c("Full", "Fusion", "Penalty", "BLasso", "BEN", "GLasso", "True")
train_beta_list_new <- lapply(methods, function(m) pad_beta_list(train_beta_list[[m]], m))
names(train_beta_list_new) <- methods

# Initialize result storage
n_covs <- length(original_levels)
metrics_names <- c("TPR", "TNR", "PPV", "NPV")

results_list <- list()

for (cov_idx in 1:n_covs) {
  # Get index combinations for this covariate
  indices <- covariate_indices[[cov_idx]]
  
  method_names <- methods[methods != "True"]
  results_cov <- matrix(NA, nrow = length(method_names), ncol = length(metrics_names))
  rownames(results_cov) <- method_names
  colnames(results_cov) <- metrics_names
  
  for (m_idx in seq_along(method_names)) {
    method <- method_names[m_idx]
    metrics_per_dataset <- matrix(NA, nrow = n_datasets, ncol = length(metrics_names))
    for (d in 1:n_datasets) {
      beta_true <- sapply(train_beta_list_new[["True"]], `[`, d)
      beta_pred <- sapply(train_beta_list_new[[method]], `[`, d)
      
      # Compute pairwise differences only for current covariate
      pairs <- combn(indices, 2)
      true_diff <- numeric(ncol(pairs))
      pred_diff <- numeric(ncol(pairs))
      for (i in 1:ncol(pairs)) {
        true_diff[i] <- as.numeric(abs(beta_true[pairs[1, i]] - beta_true[pairs[2, i]]) > threshold)
        pred_diff[i] <- as.numeric(abs(beta_pred[pairs[1, i]] - beta_pred[pairs[2, i]]) > threshold)
      }
      
      metrics_per_dataset[d, ] <- get_metrics(true_diff, pred_diff)
    }
    
    # Average across datasets
    results_cov[m_idx, ] <- round(colMeans(metrics_per_dataset, na.rm = TRUE) * 100, 1)
  }
  
  # Store results
  results_list[[cov_idx]] <- results_cov
}



final_df <- do.call(rbind, lapply(seq_along(results_list), function(i) {
  df <- as.data.frame(results_list[[i]])
  df$Method <- rownames(df)
  df$Covariate <- i
  df <- df[, c("Covariate", "Method", metrics_names)]
  return(df)
}))

# Reorder nicely
final_df <- final_df %>% arrange(Covariate, Method)

# Print the final table
print(final_df)

save(final_df,file = "table_TPR_NPR.RData")
