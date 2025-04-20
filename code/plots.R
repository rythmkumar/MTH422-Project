load("coef_mse_list.RData")
load("test_data.RData")
load("training_betas.RData")
load("training_mse.RData")
save("test_mse.RData")


# Load required libraries
library(ggplot2)
library(dplyr)

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

# For training MSE
training_data <- data.frame(Method = character(), MSE = numeric())
for (m in methods) {
  temp <- data.frame(Method = m, MSE = train_mse_list[[m]])
  rownames(temp) <- NULL
  training_data <- rbind(training_data, temp)
}

# pred_data

# Close all open graphics devices
while (!is.null(dev.list())) dev.off()


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
  ggsave(paste0("C", i, ".png"), p)
  # Optional: Pause between plots if needed
  readline(prompt = "Press [Enter] to continue to next plot...")
}


training_data$Method <- factor(training_data$Method,
                                 levels = methods)

# Boxplot for prediction MSE
plot_train_mse <- ggplot(training_data, aes(x = Method, y = MSE)) +
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

ggsave("training_mse.png", plot_train_mse)

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
